#loading pkgs
library(clusterProfiler)
library(enrichplot)
library(ggtreeExtra)
library(ggplot2)

plot_enrichment <- function(gs, title, n = 10) {
  dotplot(gs, facet = "intersect", showCategory = n,
          split = "intersect", label_format = 60) +
    scale_color_gradientn(colors = c("#e06663", "#327eba"),
      guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    guides(size = guide_legend(override.aes = list(shape = 1))) +
    theme(
      panel.grid.major.y = element_line(linetype = "dotted",
                                        color = "#808080"),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 1)
    ) +
    ggtitle(title) + xlab(NULL)
}

#IBD-only-gene
dir <- "IBD_2_subtypes_example/result"
genelist <- readRDS(file.path(dir,"IBD.gene.mpse.rds"))
gs <- compareCluster(geneClusters = genelist,
                     fun = "enrichKEGG",
                     organism = "ko")
saveRDS(gs, file = file.path(dir,"ko-ora.rds"))
p1 <- plot_enrichment(gs = gs,
  title = "Functional enrichment of intestinal genes"
)
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
p2 <- plot_enrichment(gs = gs,
  title = "Functional enrichment of chemical compounds"
)
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
