#clinical info
library(MicrobiotaProcess)
library(MicrobiotaProcess)
#2 omics clinical info
meta.mb <- read.csv("IBD_2_subtypes_example/input_data/metabolism_meta.csv",header = T)
meta.mg <- read.csv("IBD_2_subtypes_example/input_data/mg.meta.csv")

# read 2 omics expression data
metabolism <-read.csv("IBD_2_subtypes_example/input_data/metabolism_expr.csv",row.names = 1,header = 1,check.name=F)
metagenome <- read.csv("IBD_2_subtypes_example/input_data/mg.expr.csv",row.names=1,header=1,check.name=F)

#difference analysis function
DA <- function(expr,meta,group,filter.p="fdr",filter){
  mpse <- MPSE(expr)
  mpse <- mpse %>% left_join(meta,by = c('Sample' = 'sample'))
  result <- mpse %>% 
    filter(Diagnosis%in%group) %>% 
    mp_diff_analysis(
      .abundance = Abundance, 
      .group = Diagnosis,
      force = TRUE,
      relative = FALSE,
      filter.p = filter.p ) %>% 
    mp_extract_feature() %>% 
    filter(Sign_Diagnosis == filter) %>%
    pull(OTU)
}

#Metagenomic differential analysis
cd.mg <- DA(expr = metagenome, meta=meta.mg,group = c("Control","CD"),filter = "CD",filter.p = "pvalue")
uc.mg <- DA(expr = metagenome,meta = meta.mg,group = c("Control","UC"),filter = "UC",filter.p = "pvalue")
mg <- list(cd=cd.mg,uc=uc.mg)

#Save the drawing data needed
saveRDS("IBD_2_subtypes_example/result/IBD.gene.mpse.rds")


#Metabolomic differential analysis
cd.mb <- DA(expr=metabolism,meta = meta.mb ,group = c("Control","CD"),filter = "CD")
uc.mb <- DA(expr=metabolism,meta = meta.mb , group = c("Control","UC"),filter = "UC")
mb<- list(cd=cd.mb,uc=uc.mb)

#Save the drawing data needed
saveRDS("IBD_2_subtypes_example/result/IBD.gene.mpse.rds")







