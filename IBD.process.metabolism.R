#clinical info
meta <- read.csv("metabolism_meta.csv",header = T,row.names = 1)

# read metagenome data
metabolism <-read.csv("metabolism_expr.csv",row.names = 1,header = 1)

colnames(metabolism) <- sub("\\.","|",colnames(metabolism))


#The expression of the three groups was obtained
result <- vector("list", 3)
names(result) <- c("CD", "UC", "Control")
for (need in names(result)) {
  filtered <- meta[meta$Diagnosis == need, ]
  need_sample <- intersect(filtered$sample, colnames(metabolism))
  metabolism_select <- metabolism[ ,need_sample]
  result[[need]] <-  t(metabolism_select)
  
}

cd_metabolism <- result$CD
uc_metabolism <- result$UC
hc_metabolism <- result$Control

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
  results <- data.frame(name = colnames(case),pvalue = pvalue,log2FC = log2FC)
  return(results)
  
}

#UC vs Control
results <- differential_analysis(case = uc_metabolism, control = hc_metabolism)
condition <- results$pvalue < 0.05 & abs(results$log2FC) > 1
difference_metabolism_uc <- results[condition, ]


#CD vs Control
results <- differential_analysis(case = cd_metabolism, control = hc_metabolism)
condition <- results$pvalue < 0.05 & abs(results$log2FC) > 1
difference_metabolism_cd <- results[condition, ]

#Up-regulated genes in CD
cd_up <- subset(difference_metabolism_cd, log2FC > 1, select = name)
#Up-regulated genes in UC
uc_up <- subset(difference_metabolism_uc, log2FC > 1, select = name)

cd_gene_list <- cd_up$name
uc_gene_list <- uc_up$name

#The results were saved for visualization
genelist <- list(cd = cd_gene_list, uc = uc_gene_list)


saveRDS(genelist,"IBD.1.high.mb.rds")
