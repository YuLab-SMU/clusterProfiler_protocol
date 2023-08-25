#clinical info
meta <- read.csv("mg.meta.csv",header = T)

# read metagenome data
metagenome <-read.csv("mg.expr.csv",row.names = 1,header = 1)


result <- vector("list", 3)
names(result) <- c("CD", "UC", "Control")
for (need in names(result)) {
  filtered <- meta[meta$Diagnosis == need, ]
  need_sample <- intersect(filtered$sample, colnames(metagenome))
  metagenome_select <- metagenome[ ,need_sample]
  result[[need]] <-  t(metagenome_select)
  
}

cd_metagenome <- result$CD
uc_metagenome <- result$UC
hc_metagenome <- result$Control

#differential analysis
differential_analysis <- function(case, control) {
  pvalue <- rep(0, ncol(control))
  log2FC <- rep(0, ncol(control))
  for (i in 1:ncol(control)) {
    case_expr <- case[, i] %>% na.omit()
    ctrl_expr <- control[, i] %>% na.omit()
    res <- wilcox.test(case_expr, ctrl_expr)
    pvalue[i] <- res$p.value
    log2FC[i] <- log2(mean(case_expr) / mean(ctrl_expr))
  }
  results <- data.frame(name = colnames(case),pvalue = pvalue,log2FC = log2FC)
  return(results)
  
}

#UC vs Control
results <- differential_analysis(case = uc_metagenome, control = hc_metagenome)
condition <- results$pvalue < 0.05 & abs(results$log2FC) > 1
difference_metagenome_uc <- results[condition, ]


#CD vs Control
results <- differential_analysis(case = cd_metagenome, control = hc_metagenome)
condition <- results$pvalue < 0.05 & abs(results$log2FC) > 1
difference_metagenome_cd <- results[condition, ]

#Up-regulated genes in CD
cd_up <- subset(difference_metagenome_cd, log2FC > 1, select = name)
#Up-regulated genes in UC
uc_up <- subset(difference_metagenome_uc, log2FC > 1, select = name)

cd_gene_list <- cd_up$name
uc_gene_list <- uc_up$name

genelist <- list(cd = cd_gene_list, uc = uc_gene_list)

saveRDS(genelist,"IBD.1.high.mg.rds")
