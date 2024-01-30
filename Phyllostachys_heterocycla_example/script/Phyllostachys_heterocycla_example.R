# prepare db
input_dir <- "Phyllostachys_heterocycla_example/input_data"
library(tictoc)
tic("bamboo example total")
tic("prepare the expression matrix of Moso bamboo sequencing")
library(DESeq2)
library(clusterProfiler)
library(ggplot2)
library(aplot)
library(ggfun)
library(enrichplot)

counts <- read.delim(file = file.path(input_dir, "counts.txt"))

group_info <- read.delim(file = file.path(input_dir, "group_info.txt"))
toc()

tic("calculate logarithmic fold changes(log2FC) using DESeq2")

count_dds <- DESeqDataSetFromMatrix(countData = counts,
                                    colData = group_info,
                                    design = ~ group)
count_dds <- DESeq(count_dds)
toc()

tic("extract log2FC, and sort genes based on the log2FC values")
time_points <- setNames(object = c("168h", "24h", "2h"),
                        nm = c("168h_vs_0h", "24h_vs_0h", "2h_vs_0h"))
all_result <- lapply(time_points, function(time_point) {
  result <- results(count_dds, tidy = TRUE,
                    contrast = c("group", time_point, "0h"))
  setNames(object = result$log2FoldChange, nm = result$row) |>
    sort(decreasing = TRUE)
})
toc()

tic("Transcription factor enrichment analysis to identify perturbed transcription factors at different time points")

tf_db <- read.delim(
  file.path(input_dir, "annot_data/regulation_from_motif_CE_Phe.txt"),
  header = FALSE,
  colClasses = c("character", "NULL", "character", rep("NULL", 4))
) |> setNames(c("TF", "targetGene"))
perturbed_TF_result <- compareCluster(all_result, fun = "GSEA",
                                      minGSSize = 10, maxGSSize = 500,
                                      pvalueCutoff = .05,
                                      TERM2GENE = tf_db,
                                      seed = 1234)
perturbed_TF_plot <- dotplot(perturbed_TF_result,
                             showCategory = 25) +
  aes(shape = I(22)) +
  coord_flip() +
  theme_minimal() +
  theme_noxaxis() +
  xlab(NULL) +
  set_enrichplot_color(
    c("#6C8FAD", "#84ADA7", "#C7B398"),
    .fun = ggplot2::scale_fill_gradientn
  )
toc()

tic("predict the biological functions possibly regulated by the perturbed transcripiton factors")

tf_id <- unique(get_plot_data(perturbed_TF_plot, "ID")[, 1])

tf_genes <- split(tf_db$targetGene, tf_db$TF)[tf_id]
go_db <- read.delim(
  file = file.path(input_dir, "annot_data/Phe_GO_annotation.txt")
)
TF_GO_result <- compareCluster(tf_genes, fun = "enricher",
                               TERM2GENE = go_db[, c(2, 1)],
                               TERM2NAME = go_db[, c(2, 3)])

TF_GO_plot <- dotplot(TF_GO_result, by = "count",
                      showCategory = 3, label_format = 40) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1,
                                   angle = 30, size = 10)) +
  xlab(NULL) + ggtitle(NULL)

TF_GO_plot
toc()

tic("Visualization of Transcription Factor Family Annotation Information")
tf_family <- read.delim(
  file.path(input_dir, "annot_data/Phe_TF_list.txt"),
  row.names = 1
)
family_data <- subset(tf_family, Gene_ID %in% tf_id)
family_data <- family_data[order(family_data$Family), ]
family_data$Gene_ID <- factor(family_data$Gene_ID,
                              levels = family_data$Gene_ID)
cols <- 
tf_family_plot <- ggplot(data = family_data,
                        aes(x = Gene_ID, y = 1, fill = Family)) +
  geom_tile() +
  scale_fill_discrete(type=c('#B3D3AA', '#4B6B5C', '#E88A71', 
                '#DEAB76', '#CD574D', '#85C1BF', '#BF38AE', 
                '#176D84', '#7D83B7', '#4040C6', '#994B41')) +
  ggfun::theme_nothing()
toc()

tic("All in one integration to reveal transcription factor perturbation and subsequent biological effects")
fig <- insert_top(tf_family_plot, perturbed_TF_plot, height = 5) |>
  insert_bottom(TF_GO_plot, height = 50)
ggsave(fig, file = "Phyllostachys_heterocycla_example/result/figure3.png",
       width = 15, height = 10)
ggsave(fig, file = "Phyllostachys_heterocycla_example/result/figure3.pdf",
       width = 15, height = 10, dpi = "print")
toc()
toc()
