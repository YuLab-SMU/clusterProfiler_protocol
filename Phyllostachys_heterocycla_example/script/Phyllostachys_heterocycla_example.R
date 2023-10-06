library(DESeq2)
library(clusterProfiler)
library(ggplot2)
library(aplot)
library(ggfun)
library(enrichplot)


# prepare db
input_dir <- "Phyllostachys_heterocycla_example/input_data"
tf_db <- read.delim(
  file.path(input_dir, "annot_data/regulation_from_motif_CE_Phe.txt"),
  header = FALSE,
  colClasses = c("character", "NULL", "character", rep("NULL", 4))
) |> setNames(c("TF", "targetGene"))


tf_family <- read.delim(
  file.path(input_dir, "annot_data/Phe_TF_list.txt"),
  row.names = 1
)

# prepare counts
counts <- read.delim(file = file.path(input_dir, "counts.txt"))

group_info <- read.delim(file = file.path(input_dir, "group_info.txt"))

# data analysis
count_dds <- DESeqDataSetFromMatrix(countData = counts,
                                    colData = group_info,
                                    design = ~ group)
count_dds <- DESeq(count_dds)
time_points <- setNames(object = c("168h", "24h", "2h"),
                        nm = c("168h_vs_0h", "24h_vs_0h", "2h_vs_0h"))
all_result <- lapply(time_points, function(time_point) {
  result <- results(count_dds, tidy = TRUE,
                    contrast = c("group", time_point, "0h"))
  setNames(object = result$log2FoldChange, nm = result$row) |>
    sort(decreasing = TRUE)
})

# GSEA analysis
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
    c("#371ea3", "#46bac2", "#b3eebe"),
    .fun = ggplot2::scale_fill_gradientn
  )


# transcription factor function annotate
go_db <- read.delim(
  file = file.path(input_dir, "annot_data/Phe_GO_annotation.txt")
)
tf_id <- unique(get_plot_data(perturbed_TF_plot, "ID")[, 1])

tf_genes <- split(tf_db$targetGene, tf_db$TF)[tf_id]
TF_GO_result <- compareCluster(tf_genes, fun = "enricher",
                               TERM2GENE = go_db[, c(2, 1)],
                               TERM2NAME = go_db[, c(2, 3)])

TF_GO_plot <- dotplot(TF_GO_result, by = "count",
                      color = "qvalue",
                      showCategory = 3, label_format = 40) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1,
                                   angle = 30, size = 10)) +
  xlab(NULL) + ggtitle(NULL)

TF_GO_plot

# transcription factor family annotate
family_data <- subset(tf_family, Gene_ID %in% tf_id)
family_data <- family_data[order(family_data$Family), ]
family_data$Gene_ID <- factor(family_data$Gene_ID,
                              levels = family_data$Gene_ID)
tf_family_plot <- ggplot(data = family_data,
                         aes(x = Gene_ID, y = 1, fill = Family)) +
  geom_tile() +
  ggsci::scale_fill_simpsons(alpha = .6) +
  ggfun::theme_nothing()

fig <- insert_top(tf_family_plot, perturbed_TF_plot, height = 5) |>
  insert_bottom(TF_GO_plot, height = 50)
fig

ggsave(fig, file = "Phyllostachys_heterocycla_example/result/fig.png",
       width = 15, height = 10)
ggsave(fig, file = "Phyllostachys_heterocycla_example/result/fig.pdf",
       width = 15, height = 10)
