#clinical info
meta <- read.csv("mg.meta.csv", header = T)

# read metagenome data
metagenome <- read.csv("mg.expr.csv", row.names = 1, header = 1)


#get 3 groups' expression
list <- split(meta$sample, meta$Diagnosis)

cd_metagenome <- t(metagenome[, list$CD])
uc_metagenome <- t(metagenome[, list$UC])
hc_metagenome <- t(metagenome[, list$Control])


#differential analysis
differential_analysis <- function(case, control) {
  pvalue <- rep(0, ncol(control))
  log2FC <- rep(0, ncol(control))
  for (i in 1:ncol(control)) {
    case_expr <- na.omit(case[, i])
    ctrl_expr <- na.omit(control[, i])
    res <- wilcox.test(case_expr, ctrl_expr)
    pvalue[i] <- res$p.value
    log2FC[i] <- log2(mean(case_expr) / mean(ctrl_expr))
  }
  results <-
    data.frame(name = colnames(case),
               pvalue = pvalue,
               log2FC = log2FC)
  return(results)
  
}

#UC vs Control
results <-differential_analysis(case = uc_metagenome, control = hc_metagenome)
condition <- results$pvalue < 0.05 & abs(results$log2FC) > 1
difference_metagenome_uc <- results[condition,]


#CD vs Control
results <-
  differential_analysis(case = cd_metagenome, control = hc_metagenome)
condition <- results$pvalue < 0.05 & abs(results$log2FC) > 1
difference_metagenome_cd <- results[condition,]

#Up-regulated genes in CD
cd_up <- subset(difference_metagenome_cd, log2FC > 1, select = name)
#Up-regulated genes in UC
uc_up <- subset(difference_metagenome_uc, log2FC > 1, select = name)

cd_gene_list <- cd_up$name
uc_gene_list <- uc_up$name

genelist <- list(cd = cd_gene_list, uc = uc_gene_list)

saveRDS(genelist, "IBD.1.high.mg.rds")
