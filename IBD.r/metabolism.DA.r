setwd("~/result/real-end/")
library(MicrobiotaProcess)
#clinical info
meta <- read.csv("metabolism_meta.csv",header = T)

# read metabolism  data
metabolism <-read.csv("metabolism_expr.csv",row.names = 1,header = 1)

colnames(metabolism) <- sub("\\.","|",colnames(metabolism))

mpse <- MPSE(metabolism )
mpse <- mpse %>% dplyr::left_join(meta, by = c('Sample' = 'sample'))

mpse.cd <- mpse %>% dplyr::filter(Diagnosis %in% c('Control', 'CD'))
mpse.uc <- mpse %>% dplyr::filter(Diagnosis %in% c('Control', 'UC'))

mpse.cd <- mpse.cd %>% 
  mp_diff_analysis(
    .abundance = Abundance, 
    .group = Diagnosis,
    force = TRUE,
    relative = FALSE
  )

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

saveRDS(genelist, "IBD.cpd.mpse.rds")
