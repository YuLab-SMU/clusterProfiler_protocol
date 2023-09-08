#load pkgs
library(clusterProfiler)
library(enrichplot)
library(ggtreeExtra)
library(ggplot2)
plot.enrichment <- function(gs, title, n) {
  dotplot(
    gs,
    facet = "intersect",
    showCategory = n,
    split = 'intersect',
    label_format = 60
  ) + scale_color_gradientn(
    colours = c("#b3eebe", "#46bac2", "#371ea3"),
    guide = guide_colorbar(reverse = TRUE, order = 1)
  ) +
    guides(size = guide_legend(override.aes = list(shape = 1))) +
    theme(
      panel.grid.major.y = element_line(linetype = 'dotted', color = '#808080'),
      panel.grid.major.x = element_blank()
    ) + ggtitle(title)
}

#IBD-only-gene
genelist<- readRDS("IBD_2_subtypes_example/result/IBD.gene.mpse.rds")
gs <- compareCluster(geneClusters = genelist , fun = "enrichKEGG", 
                     organism = "ko")
plot.enrichment(gs = gs , title = "IBD 2 subtypes gene ORA", n = 15)
ggsave("IBD_2_subtypes_example/result/IBD_2_subtypes_gene_ORA.pdf",
       width = 9,
       height = 9)

gs2  <-  pairwise_termsim(gs)
treeplot(gs2, offset_tiplab = 9) + ggtitle("IBD 2 subtypes gene ORA")
ggsave(
  "IBD_2_subtypes_example/result/IBD_2_subtypes_in_metagenome_ORA(treeplot).pdf",
  width = 15,
  height = 5.5
)


#IBD-only-compound
cpd.list <-  readRDS("IBD_2_subtypes_example/result/IBD.cpd.mpse.rds")
gs <- compareCluster(geneClusters = cpd.list, fun = "enrichKEGG", organism = "cpd")
plot.enrichment(gs = gs, title = "IBD 2 subtypes compound ORA", n = 15)
ggsave("IBD_2_subtypes_example/result/IBD_2_subtypes_in_metabolism_ORA.pdf", width = 9, height = 9)

gs2  <-  pairwise_termsim(gs)
treeplot(gs2, offset_tiplab = 10, showCategory = 15) + ggtitle("IBD 2 subtypes compound ORA")
ggsave("IBD_2_subtypes_example/result/IBD_2_subtypes_in_metabolism_ORA(treeplot).pdf", width = 17, height = 12)


