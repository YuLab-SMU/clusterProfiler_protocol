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
                       check.name = FALSE)

metagenome <- read.csv(file.path(input_dir, "mg.expr.csv"),
                       row.names = 1,
                       header = TRUE,
                       check.name = FALSE)

# difference analysis function
DA <- function(expr, 
               meta, 
               abundance = 'Abundance', 
               group.colname = 'Diagnosis', 
               force = TRUE,
               relative = FALSE,
               subset.group, 
               diff.group, 
               filter.p = 'pvalue', ...) {
  sign.group.colname <- paste0('Sign_', group.colname)
  mpse <- MPSE(expr)
  mpse <- mpse |> left_join(meta, by = "Sample")
  mpse |>
    dplyr::filter(!!as.symbol(group.colname) %in% subset.group) |>
    mp_diff_analysis(
      .abundance = !!as.symbol(abundance),
      .group = !!as.symbol(group.colname),
      force = force,
      relative = relative,
      filter.p = filter.p,
      ...
    ) |>
    mp_extract_feature() |>
    dplyr::filter(!!as.symbol(sign.group.colname) == diff.group) |>
    dplyr::pull(OTU)
}

# Metagenomic differential analysis
groups <- c(CD = "CD", UC = "UC")
to.ora.gene <- lapply(groups, function(x){
  DA(expr = metagenome,
    meta = meta_mg,
    subset.group = c(x, 'Control'), 
    diff.group = x,
  )
})

# Save the drawing data needed
saveRDS(to.ora.gene, file.path(output_dir, "IBD.gene.mpse.rds"))

# Metabolomic differential analysis
groups <- c(CD = "CD", UC = "UC")
to.ora.cpd <- lapply(groups, function(x){
  DA(expr = metabolism,
    meta = meta_mb,
    subset.group = c(x, "Control"),
    diff.group = x,
  )
})

# Save the drawing data needed
saveRDS(to.ora.cpd, file.path(output_dir, "IBD.cpd.mpse.rds"))
