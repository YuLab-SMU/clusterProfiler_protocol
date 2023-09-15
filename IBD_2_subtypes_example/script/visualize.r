#loading pkgs
library(clusterProfiler)
library(enrichplot)
library(ggtreeExtra)
library(ggplot2)




plot_enrichment <- function(gs, title, n = 10) {
  dotplot(gs, facet = "intersect", showCategory = n,
          split = "intersect", label_format = 60) +
    theme(
      panel.grid.major.y = element_line(linetype = "dotted",
                                        color = "#808080"),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 1)
    ) +
    np_style()
}

#IBD-only-gene
dir <- "IBD_2_subtypes_example/result"
genelist <- readRDS(file.path(dir,"IBD.gene.mpse.rds"))
gs <- compareCluster(geneClusters = genelist,
                     fun = "enrichKEGG",
                     organism = "ko")
saveRDS(gs, file = file.path(dir,"ko-ora.rds"))

gs <- readRDS(file.path(dir,"ko-ora.rds"))

p1 <- plot_enrichment(gs = gs,
  title = "Functional enrichment of intestinal genes"
)

p1 <- dotplot(gs, facet='intersect', showCategory = 10, split = "intersect", label_format = 60) +
  ggtitle("Functional enrichment of intestinal genes") + 
  theme(plot.title = element_text(hjust = 1))

ggsave(p1,
       filename = file.path(dir, "IBD_2_subtypes_gene_ORA.pdf"),
       width = 9,
       height = 9)

#IBD-only-compound
cpd_list <- readRDS(file.path(dir, "IBD.cpd.mpse.rds"))
gs <- compareCluster(geneClusters = cpd_list,
                     fun = "enrichKEGG",
                     organism = "cpd")
saveRDS(gs, file = file.path(dir, "cpd-ora.rds"))
gs <- readRDS(file.path(dir, "cpd-ora.rds"))


p2 <- plot_enrichment(gs = gs,
  title = "Functional enrichment of chemical compounds"
)

p2 <- dotplot(gs, facet='intersect', showCategory = 10, split = "intersect", label_format = 60) +
  ggtitle("Functional enrichment of chemical compounds") +
  theme(plot.title = element_text(hjust = 1))

ggsave(p2,
       filename = file.path(dir,"IBD_2_subtypes_in_metabolism_ORA.pdf"),
       width = 9,
       height = 9)

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
