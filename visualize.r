## install.packages
devtools::install_github("YuLab-SMU/enrichplot")

#load pkgs
library(clusterProfiler)
library(gson)
library(enrichplot)
library(ggplot2)
library(ggtreeExtra)

#two gson

gson_cpd <- function(){
  k1 <-clusterProfiler::: kegg_rest("https://rest.kegg.jp/link/cpd/pathway")
  k1[,1] %<>% gsub("[^:]+:", "", .)
  k1[,2] %<>% gsub("[^:]+:", "", .)
  k1 <- k1[grep("map",k1[,1]),]
  k2  <- clusterProfiler:::kegg_rest("https://rest.kegg.jp/list/pathway")
  k2[,1] %<>% gsub("path:","",.)
  gsid2gene <- setNames(k1, c("gsid", "gene"))
  gsid2name <- setNames(k2, c("gsid", "name"))
  y <- readLines("https://rest.kegg.jp/info/ko")
  version <- sub("\\w+\\s+", "", y[grep('Release', y)])
  gson::gson(gsid2gene = gsid2gene,
             gsid2name = gsid2name,
             species = "KEGG Compound",
             gsname = "KEGG",
             version = version,
             keytype = "kegg_compound",
             accessed_date = as.character(Sys.Date()))
}

gson_KO <- function() {
    k1 <- clusterProfiler:::kegg_rest("https://rest.kegg.jp/link/ko/pathway")
    k1[,1] %<>% gsub("[^:]+:", "", .)
    k1[,2] %<>% gsub("[^:]+:", "", .)
    k1 <- k1[grep("map",k1[,1]),]
    k2  <- clusterProfiler:::kegg_rest("https://rest.kegg.jp/list/pathway")
    k2[,1] %<>% gsub("path:","",.)
    gsid2gene <- setNames(k1, c("gsid", "gene"))
    gsid2name <- setNames(k2, c("gsid", "name"))
    y <- readLines("https://rest.kegg.jp/info/ko")
    version <- sub("\\w+\\s+", "", y[grep('Release', y)])
    gson(gsid2gene = gsid2gene,
         gsid2name = gsid2name,
         species = "KEGG Orthology",
         gsname = "KEGG",
         version = version,
         keytype = "kegg_orthology",
         accessed_date = as.character(Sys.Date()))
}
plot.enrichment <- function(gs,title,n){
dotplot(gs, facet = "intersect", showCategory=n,split='set')+
 scale_color_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
 guide=guide_colorbar(reverse=TRUE, order=1)) +
 guides(size = guide_legend(override.aes=list(shape=1))) +
 theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
 panel.grid.major.x = element_blank())+ggtitle(title)
}

#IBD only
#IBD-only-gene
genelist <- readRDS("IBD.1.high.mg.rds")
gs=compareCluster(geneClusters = genelist,fun="enricher",gson = gson_KO())

#dotplot
plot.enrichment(gs=gs ,title = "IBD 2 subtype gene ORA",n=15)
#tree plot
gs2=pairwise_termsim(gs)
treeplot(gs2,offset_tiplab=7,showCategory=15)+ggtitle("IBD 2 subtype gene ORA")

#IBD-only-compound
cpd.list <- readRDS("IBD.1.high.mb.rds")
gs=compareCluster(geneClusters = cpd.list,fun="enricher",gson = gson_cpd())

#dotplot
plot.enrichment(gs=gs,title = "IBD 2 subtype compound ORA",n=15)

#tree plot
gs2=pairwise_termsim(gs)
treeplot(gs2,offset_tiplab = 5,showCategory=15)+ggtitle("IBD 2 subtype compound ORA")





