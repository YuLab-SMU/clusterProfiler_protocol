#loading pkgs
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
# devtools::install_github("YuLab-SMU/enrichplot")
#IBD-only-gene
dir <- "IBD_2_subtypes_example/result"
de_gene <- readRDS(file.path(dir,"IBD.gene.mpse.rds"))
gs <- compareCluster(geneClusters = de_gene,
                     fun = "enrichKEGG",
                     organism = "ko")


p1 <- dotplot(gs, facet='intersect', showCategory = 10, split = "intersect", label_format = 60) +
  ggtitle("Functional enrichment of intestinal genes") + 
  theme(plot.title = element_text(hjust = 1))

ggsave(p1,
       filename = file.path(dir, "IBD_2_subtypes_gene_ORA.pdf"),
       width = 9,
       height = 9)

#IBD-only-compound
de_cpd <- readRDS(file.path(dir, "IBD.cpd.mpse.rds"))
gs <- compareCluster(geneClusters = de_cpd,
                     fun = "enrichKEGG",
                     organism = "cpd")


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
