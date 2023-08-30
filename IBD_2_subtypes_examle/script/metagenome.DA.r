library(MicrobiotaProcess)

#clinical info
meta <- read.csv("../input_data/mg.meta.csv", header = T)

# read metagenome data
metagenome <- read.csv("../input_data/mg.expr.csv", row.names = 1, header = 1)
colnames(metagenome) <- sub("\\.", "|", colnames(metagenome))
mpse <- MPSE(metagenome)

#clinical info
mpse <- mpse %>% dplyr::left_join(meta, by = c('Sample' = 'sample'))
mpse.cd <- mpse %>% dplyr::filter(Diagnosis %in% c('Control', 'CD'))
mpse.uc <- mpse %>% dplyr::filter(Diagnosis %in% c('Control', 'UC'))

#CD vs Control
mpse.cd.res<- mpse.cd %>% 
  mp_diff_analysis(
    .abundance = Abundance, 
    .group = Diagnosis,
    force = TRUE,
    relative = FALSE,
    filter.p = "pvalue"
    
  )
#UC vs Control
mpse.uc.res <- mpse.uc %>%
  mp_diff_analysis(
    .abundance = Abundance,
    .group = Diagnosis,
    force = TRUE,
    relative = FALSE,
    filter.p = "pvalue"
    
  )




genelist <- list(cd = mpse.cd.res %>% 
                   mp_extract_feature() %>% 
                   dplyr::filter(Sign_Diagnosis == 'CD') %>%
                   dplyr::pull(OTU),

                 uc = mpse.uc.res %>%
                   mp_extract_feature() %>%
                   dplyr::filter(Sign_Diagnosis == 'UC') %>%
                   dplyr::pull(OTU)
)

saveRDS(genelist, "../result/IBD.gene.mpse.rds")
