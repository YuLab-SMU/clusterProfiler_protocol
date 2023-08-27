#load pkgs
library(clusterProfiler)
library(gson)
library(enrichplot)
library(ggplot2)
library(ggtreeExtra)

#two gson
kegg_rest <- getFromNamespace("kegg_rest", "clusterProfiler")

gson_cpd <- function() {
  k1 <- kegg_rest("https://rest.kegg.jp/link/cpd/pathway")
  k1[, 1]  <- gsub("[^:]+:", "", k1[, 1])
  k1[, 2]  <-  gsub("[^:]+:", "",  k1[, 2])
  k1 <- k1[grep("map", k1[, 1]),]
  k2  <- kegg_rest("https://rest.kegg.jp/list/pathway")
  k2[, 1]  <- gsub("path:", "", k2[, 1])
  gsid2gene <- setNames(k1, c("gsid", "gene"))
  gsid2name <- setNames(k2, c("gsid", "name"))
  y <- readLines("https://rest.kegg.jp/info/ko")
  version <- sub("\\w+\\s+", "", y[grep('Release', y)])
  gson::gson(
    gsid2gene = gsid2gene,
    gsid2name = gsid2name,
    species = "KEGG Compound",
    gsname = "KEGG",
    version = version,
    keytype = "kegg_compound",
    accessed_date = as.character(Sys.Date())
  )
}

gson_KO <- function() {
  k1 <- kegg_rest("https://rest.kegg.jp/link/ko/pathway")
  k1[, 1]  <- gsub("[^:]+:", "", k1[, 1])
  k1[, 2]  <-  gsub("[^:]+:", "",  k1[, 2])
  k1 <- k1[grep("map", k1[, 1]), ]
  k2  <- kegg_rest("https://rest.kegg.jp/list/pathway")
  k2[, 1]  <- gsub("path:", "", k2[, 1])
  gsid2gene <- setNames(k1, c("gsid", "gene"))
  gsid2name <- setNames(k2, c("gsid", "name"))
  y <- readLines("https://rest.kegg.jp/info/ko")
  version <- sub("\\w+\\s+", "", y[grep('Release', y)])
  gson(
    gsid2gene = gsid2gene,
    gsid2name = gsid2name,
    species = "KEGG Orthology",
    gsname = "KEGG",
    version = version,
    keytype = "kegg_orthology",
    accessed_date = as.character(Sys.Date())
  )
}
plot.enrichment <- function(gs, title, n) {
  dotplot(gs,facet = "intersect",
          showCategory = n, split = 'set') +
    scale_color_gradientn( colours = c("#b3eebe", "#46bac2", "#371ea3"),
                           guide = guide_colorbar(reverse = TRUE, order =1)) +
    guides(size = guide_legend(override.aes = list(shape = 1))) +
    theme(panel.grid.major.y = element_line(linetype = 'dotted', color = '#808080'),
          panel.grid.major.x = element_blank()) + ggtitle(title)
}

#IBD-only-gene
genelist<- readRDS("IBD.gene.mpse.rds")
gs <- compareCluster(geneClusters = genelist , fun = "enricher",gson = gson_KO())
plot.enrichment(gs = gs , title = "IBD 2 subtypes gene ORA", n = 15)
ggsave("IBD_2_subtypes_gene_ORA.pdf",width =5.5,height=7 )

gs2 = pairwise_termsim(gs)
treeplot(gs2, offset_tiplab =9) + ggtitle("IBD 2 subtypes gene ORA")
ggsave("IBD_2_subtypes_in_metagenome_ORA(treeplot).pdf", width = 15, height = 5.5)


#IBD-only-compound
cpd.list <-  readRDS("IBD.cpd.mpse.rds")
gs = compareCluster(geneClusters = cpd.list, fun = "enricher", gson = gson_cpd())
plot.enrichment(gs = gs, title = "IBD 2 subtypes compound ORA", n = 15)
ggsave("IBD_2_subtypes_in_metabolism_ORA.pdf", width = 5.5, height = 7)
enricher()

gs2 = pairwise_termsim(gs)
treeplot(gs2, offset_tiplab = 10, showCategory = 15) + ggtitle("IBD 2 subtypes compound ORA")
ggsave("IBD_2_subtypes_in_metabolism_ORA(treeplot).pdf", width = 40, height = 5.5)
