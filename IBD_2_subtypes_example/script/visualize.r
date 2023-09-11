#loading pkgs
library(clusterProfiler)
library(enrichplot)
library(ggtreeExtra)
library(ggplot2)

plot.enrichment <- function(gs, title, n = 10) {
  dotplot(
    gs,
    facet = "intersect",
    showCategory = n,
    split = 'intersect',
    label_format = 60
  ) +
    scale_color_gradientn(
      colors = c("#e06663", "#327eba"),
      guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    guides(size = guide_legend(override.aes = list(shape = 1))) +
    theme(
      panel.grid.major.y = element_line(linetype = 'dotted', color = '#808080'),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 1)
    ) +
    ggtitle(title) + xlab(NULL)
}

#IBD-only-gene
genelist <- readRDS("IBD_2_subtypes_example/result/IBD.gene.mpse.rds")
gs <- compareCluster(geneClusters = genelist ,
                     fun = "enrichKEGG",
                     organism = "ko")
saveRDS(gs, file = "IBD_2_subtypes_example/result/ko-ora.rds")
p1 <- plot.enrichment(gs = gs , title = "Functional enrichment of intestinal genes")
ggsave(p1,
       filename = "IBD_2_subtypes_example/result/IBD_2_subtypes_gene_ORA.pdf",
       width = 9,
       height = 9)

#IBD-only-compound
cpd.list <- readRDS("IBD_2_subtypes_example/result/IBD.cpd.mpse.rds")
gs <- compareCluster(geneClusters = cpd.list,
                 fun = "enrichKEGG",
                 organism = "cpd")
saveRDS(gs, file = "IBD_2_subtypes_example/result/cpd-ora.rds")
p2 <- plot.enrichment(gs = gs, title = "Functional enrichment of chemical compounds")
ggsave(p2,
       filename = "IBD_2_subtypes_example/result/IBD_2_subtypes_in_metabolism_ORA.pdf",
       width = 9,
       height = 9)

fig = aplot::plot_list(p1, p2, tag_levels = 'A', tag_size = 15,
                       widths = c(1, .8))
ggsave(fig,
       file = "IBD_2_subtypes_example/result/fig.pdf",
       width = 13,
       height = 6.5)
ggsave(fig,
       file = "IBD_2_subtypes_example/result/fig.png",
       width = 13,
       height = 6.5)