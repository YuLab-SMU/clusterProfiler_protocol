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
DA <- function(expr, meta, abundance, group.colname, subset.group, 
               diff.group, filter.p = 'pvalue', ...) {
  sign.group.colname <- paste0('Sign_', group.colname)
  mpse <- MPSE(expr)
  mpse <- mpse |> left_join(meta, by = "Sample")
  mpse |>
    dplyr::filter(!!as.symbol(group.colname) %in% subset.group) |>
    mp_diff_analysis(
      .abundance = !!as.symbol(abundance),
      .group = !!as.symbol(group.colname),
      force = TRUE,
      relative = FALSE,
      filter.p = filter.p,
      ...
    ) |>
    mp_extract_feature() |>
    dplyr::filter(!!as.symbol(sign.group.colname) == diff.group) |>
    dplyr::pull(OTU)
}

# Metagenomic differential analysis
subset.groups <- list(c('Control', 'CD'), c('Control', 'UC'))
mg <- lapply(subset.groups, function(x){
    DA(expr = metagenome, 
       meta = meta_mg, 
       abundance = "Abundance",
       group.colname = "Diagnosis", 
       subset.group = x, 
       diff.group = setdiff(x, 'Control'),
    )
})
names(mg) <- c('CD', 'UC')

# Save the drawing data needed
saveRDS(mg, file.path(output_dir, "IBD.gene.mpse.rds"))

# Metabolomic differential analysis
mb <- lapply(subset.groups, function(x){
    DA(expr = metabolism,
       meta = meta_mb,
       abundance = "Abundance",
       group.colname = 'Diagnosis',
       subset.group = x,
       diff.group = setdiff(x, 'Control'),
    )
})
names(mb) <- c('CD', 'UC')

# Save the drawing data needed
saveRDS(mb, file.path(output_dir, "IBD.cpd.mpse.rds"))
