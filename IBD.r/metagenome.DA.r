setwd("~/result/real-end/")
library(MicrobiotaProcess)

#clinical info
meta <- read.csv("mg.meta.csv", header = T)

# read metagenome data
metagenome <- read.csv("mg.expr.csv", row.names = 1, header = 1)
colnames(metagenome) <- sub("\\.", "|", colnames(metagenome))
mpse <- MPSE(metagenome)

#clinical info
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


# UC vs Control
p3 <-
  mpse.uc %>% as_tibble() %>% as_tibble()  %>% dplyr::group_by(OTU)  %>% tidytable::group_split(OTU, .named =
                                                                                                  T)  %>%
  purrr::map(function(x) {
    wilcox.test(Abundance ~ Diagnosis, data = as.data.frame(x))$p.value
  }) %>% unlist()

names <- p3[p3 < 0.05] %>% na.omit()  %>% names()
mpse.uc2 <- mpse.uc %>% filter(OTU %in% names)
mpse.uc3 <- mpse.uc2 %>% dplyr::group_by(OTU, Diagnosis) %>%
  dplyr::summarize(mean.abundance = mean(Abundance + 1)) %>%
  tidytable::group_split(OTU, .named = T) %>% purrr::map(function(x)
    log2(x[2, 3] / x[1, 3])) %>%
  dplyr::bind_rows(.id = 'OTU') %>% dplyr::rename(log2fc = mean.abundance) %>% filter(log2fc >
                                                                                        1)


genelist <- list(
  cd = mpse.cd %>%
    mp_extract_feature() %>%
    dplyr::filter(Sign_Diagnosis == 'CD') %>%
    dplyr::pull(OTU),
  
  uc = mpse.uc3$OTU
)

saveRDS(genelist, "IBD.gene.mpse.rds")
