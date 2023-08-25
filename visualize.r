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
  k1 <- k1[grep("map", k1[, 1])]
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
#IBD onl
#IBD-only-gene
genelist <- readRDS("~/result/nature-protocal-818/new-code/IBD/IBD.1.high.mg.rds")
gs = compareCluster(geneClusters = genelist, fun = "enricher",gson = gson_KO())
pdf("IBD_2_subtype_in_metagenome_ORA.pdf", width = 5.5, height = 7)
plot.enrichment(gs = gs , title = "IBD 2 subtype gene ORA", n = 15)
dev.off()

pdf("IBD_2_subtype_in_metagenome_ORA(treeplot).pdf", width = 13, height = 5.5)
gs2 = pairwise_termsim(gs)
treeplot(gs2, offset_tiplab = 10) + ggtitle("IBD 2 subtype gene ORA")
dev.off()

#IBD-only-compound
cpd.list <- readRDS("~/result/nature-protocal-818/new-code/IBD/IBD.1.high.mb.rds")
gs = compareCluster(geneClusters = cpd.list, fun = "enricher", gson = gson_cpd())

pdf("IBD_2_subtype_in_metabolism_ORA.pdf", width = 5.5, height = 7)
plot.enrichment(gs = gs, title = "IBD 2 subtype compound ORA", n = 15)
dev.off()


pdf("IBD_2_subtype_in_metabolism_ORA(treeplot).pdf", width = 13, height = 5.5)
gs2 = pairwise_termsim(gs)
treeplot(gs2, offset_tiplab = 5, showCategory = 15) + ggtitle("IBD 2 subtype compound ORA")
dev.off()
