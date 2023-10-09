library(tictoc)
dir <- "IBD_2_subtypes_example/result"
tic("Procedure 1: Metabolomics and Metagenomics Functional Enrichment Analysis")
tic("setup the environment and data objects")
library(MicrobiotaProcess)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)

meta_mg <- read.csv("IBD_2_subtypes_example/input_data/mg.meta.csv")
metagenome <- read.csv("IBD_2_subtypes_example/input_data/mg.expr.csv",
                       row.names = 1,
                       check.name = FALSE)
toc()


tic("Metagenomic Data Differential Analysis")
DA <- function(expr,
               meta,
               abundance = 'Abundance',
               group_colname = 'Diagnosis',
               force = TRUE,
               relative = FALSE,
               subset_group,
               diff_group,
               filter.p = 'pvalue', ...) {
  sign_group_colname <- paste0('Sign_', group_colname)
  mpse <- MPSE(expr)
  mpse <- mpse |> left_join(meta, by = "Sample")
  mpse |>
    dplyr::filter(!!as.symbol(group_colname) %in% subset_group) |>
    mp_diff_analysis(
      .abundance = !!as.symbol(abundance),
      .group = !!as.symbol(group_colname),
      force = force,
      relative = relative,
      filter.p = filter.p,
      ...
    ) |>
    mp_extract_feature() |>
    dplyr::filter(!!as.symbol(sign_group_colname) == diff_group) |>
    dplyr::pull(OTU)
}

groups <- c(CD = 'CD', UC = 'UC')
de_gene <- lapply(groups, function(x) {
  DA(expr = metagenome,
     meta = meta_mg,
     subset_group = c(x, 'Control'),
     diff_group = x)
})
toc()

tic("Functional analysis of differential microbial genes")
gene_enrich_result <- compareCluster(geneClusters = de_gene,
                                     fun = "enrichKEGG",
                                     organism = "ko")

p1 <- dotplot(gene_enrich_result, facet = 'intersect', showCategory = 10,
        split = "intersect", label_format = 60) +
  ggtitle("Functional enrichment of intestinal genes") +
  theme(plot.title = element_text(hjust = 1))
ggsave(p1,
       filename = file.path(dir, "IBD_2_subtypes_gene_ORA.pdf"),
       width = 9,
       height = 9)
toc()

tic("Metabolomic data differential analysis")
meta_mb <- read.csv("IBD_2_subtypes_example/input_data/metabolism_meta.csv")
metabolism <- read.csv("IBD_2_subtypes_example/input_data/metabolism_expr.csv",
                       row.names = 1,
                       check.name = FALSE)

groups <- c(CD = 'CD', UC = 'UC')
de_cpd <- lapply(groups, function(x) {
  DA(expr = metabolism,
     meta = meta_mb,
     subset_group = c(x, 'Control'),
     diff_group = x)
})
toc()

tic("Functional analysis of differential metabolites")
cpd_enrich_result <- compareCluster(geneClusters = de_cpd,
                                    fun = "enrichKEGG",
                                    organism = "cpd")

p2 <- dotplot(cpd_enrich_result, facet = 'intersect', showCategory = 10,
        split = "intersect", label_format = 60) +
  ggtitle("Functional enrichment of chemical compounds") +
  theme(plot.title = element_text(hjust = 1))
ggsave(p2,
       filename = file.path(dir,"IBD_2_subtypes_in_metabolism_ORA.pdf"),
       width = 9,
       height = 9)
toc()
toc()

fig  <-  aplot::plot_list(p1, p2, tag_levels = "A", tag_size = 15,
                          widths = c(1, .8))
ggsave(fig,
       file = file.path(dir, "fig.pdf"),
       width = 13,
       height = 6.5)
ggsave(fig,
       file = file.path(dir, "fig.png"),
       width = 13,
       height = 6.5)