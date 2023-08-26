library(MicrobiotaProcess)

#clinical info
meta <- read.csv("mg.meta.csv", header = T)

# read metagenome data
metagenome <- read.csv("mg.expr.csv", row.names = 1, header = 1)

mpse <- MPSE(metagenome)
mpse <- mpse %>% dplyr::left_join(meta, by = c('Sample' = 'sample'))

mpse.cd <- mpse %>% dplyr::filter(Diagnosis %in% c('Control', 'CD'))
mpse.uc <- mpse %>% dplyr::filter(Diagnosis %in% c('Control', 'UC'))

#CD vs Control
mpse.cd <- mpse.cd %>% 
  mp_diff_analysis(
    .abundance = Abundance, 
    .group = Diagnosis,
    force = TRUE,
    relative = FALSE
  )
#UC vs Control
mpse.uc <- mpse.uc %>%
  mp_diff_analysis(
    .abundance = Abundance,
    .group = Diagnosis,
    force = TRUE,
    relative = FALSE
  )


genelist <- list(cd = mpse.cd %>% 
                   mp_extract_feature() %>% 
                   dplyr::filter(Sign_Diagnosis == 'CD') %>%
                   dplyr::pull(OTU),
                 
                 uc = mpse.uc %>%
                   mp_extract_feature() %>%
                   dplyr::filter(Sign_Diagnosis == 'UC') %>%
                   dplyr::pull(OTU)
)

saveRDS(genelist, "IBD.gene.mpse.rds")
