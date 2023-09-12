# loading pkg
library(MicrobiotaProcess)
# 2 omics clinical info
input_dir <- "IBD_2_subtypes_example/input_data"
output_dir <- "IBD_2_subtypes_example/result"
meta_mb <- read.csv(file.path(input_dir, "metabolism_meta.csv"), header = TRUE)
meta_mg <- read.csv(file.path(input_dir, "mg.meta.csv"), header = TRUE)

# read 2 omics expression data
metabolism <- read.csv(file.path(input_dir, "metabolism_expr.csv"),
                       row.names = 1,
                       header = TRUE,
                       check.name = FALSE
)
metagenome <- read.csv(file.path(input_dir, "mg.expr.csv"),
                       row.names = 1,
                       header = TRUE,
                       check.name = FALSE
)

# difference analysis function
DA <- function(expr, meta, diff_group) {
  mpse <- MPSE(expr)
  mpse <- mpse %>% left_join(meta, by = "Sample")
  mpse %>%
    filter(Diagnosis %in% c("Control", diff_group)) %>%
    mp_diff_analysis(
      .abundance = Abundance,
      .group = Diagnosis,
      force = TRUE,
      relative = FALSE,
      filter.p = "pvalue"
    ) %>%
    mp_extract_feature() %>%
    filter(Sign_Diagnosis == diff_group) %>%
    pull(OTU)
}

# Metagenomic differential analysis
cases <- c("CD", "UC")
mg <- lapply(cases, function(case) {
  DA(expr = metagenome, meta = meta_mg, diff_group = case)
})
names(mg) <- cases
# Save the drawing data needed
saveRDS(mg, file.path(output_dir, "IBD.gene.mpse.rds"))

# Metabolomic differential analysis
mb <- lapply(cases, function(case) {
  DA(expr = metabolism, meta = meta_mb, diff_group = case)
})
names(mb) <- cases
# Save the drawing data needed
saveRDS(mb, file.path(output_dir, "IBD.cpd.mpse.rds"))