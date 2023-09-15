library(DESeq2)
library(clusterProfiler)
library(ggplot2)
library(aplot)
library(ggfun)
library(enrichplot)

#np_style <- function(colors = c("#e06663", "#327eba"), legend_shape=1, title = NULL) {
#  list(scale_color_gradientn(colors = colors,
#        guide = guide_colorbar(reverse = TRUE, order = 1)),
#    guides(size = guide_legend(override.aes = list(shape = legend_shape))),
#    ggtitle(title),
#    xlab(NULL)
#  )
#}


# prepare db
input_dir <- "Phyllostachys_heterocycla_example/input_data"
plant_tf_db <- read.table(
  file.path(input_dir, "annot_data/regulation_from_motif_CE_Phe.txt"),
  sep = "\t", header = FALSE
)
plant_tf_db <- plant_tf_db[, c(1, 3)]
colnames(plant_tf_db) <- c("TERM", "GENE")

tf_id_annotation <- read.table(
  file.path(input_dir, "annot_data/Phe_TF_list.txt"),
  header = TRUE, sep = "\t"
)

# prepare counts
counts <- read.table(
  file = file.path(input_dir, "counts.txt"),
  header = TRUE, sep = "\t"
)

group_info <- read.table(
  header = TRUE,
  file = file.path(input_dir, "group_info.txt")
)

# data analysis
count_dds <- DESeqDataSetFromMatrix(countData = counts,
                                    colData = group_info,
                                    design = ~group)
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
set.seed(1234)
compare_enrich_result <- compareCluster(all_result, fun = "GSEA",
                                        minGSSize = 10, maxGSSize = 500,
                                        pvalueCutoff = .05,
                                        TERM2GENE = plant_tf_db)

if (FALSE) {
enrich_plot <- dotplot(compare_enrich_result, includeAll = TRUE,
                       showCategory = 25) +
  aes(shape = I(15)) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(linetype = "dotted",
                                          color = "#808080"),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank()) +
  np_style(colours = c("#371ea3", "#46bac2", "#b3eebe"), legend_shape=0)
}



enrich_plot <- dotplot(compare_enrich_result, showCategory = 25) + 
  aes(shape=I(22)) +
  coord_flip() +
  theme_minimal() +
  theme_noxaxis() + 
  xlab(NULL) +
  set_enrichplot_color(c("#371ea3", "#46bac2", "#b3eebe"), .fun=ggplot2::scale_fill_gradientn)



# transcription factor function annotate
go_db <- read.table(
  file = file.path(input_dir, "annot_data/Phe_GO_annotation.txt"),
  header = TRUE, sep = "\t", quote = ""
)

# tf_id <- unique(enrich_plot$data$ID)
tf_id <- unique(get_plot_data(enrich_plot, "ID")[,1])

tf_genes <- split(plant_tf_db$GENE, plant_tf_db$TERM)[tf_id]
y <- compareCluster(tf_genes, fun = "enricher",
                    TERM2GENE = go_db[, c(2, 1)],
                    TERM2NAME = go_db[, c(2, 3)])

if (FALSE) {
enrich_pathway_plot <- dotplot(y, by = "count", color = "qvalue",
                               showCategory = 3, label_format = 40) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1,
                                   angle = 30, size = 10),
        axis.text.y = element_text(size = 13, face = "bold")) +
  np_style()
}

enrich_pathway_plot <- dotplot(y, by = "count", color = "qvalue",
          showCategory = 3, label_format = 40) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1,
                                   angle = 30, size = 10)) +
  xlab(NULL) + ggtitle(NULL)

enrich_pathway_plot
#  guides(size = guide_legend(override.aes = list(shape = 1))) +
#  xlab(NULL) + ggtitle(NULL) +
#  scale_color_gradientn(colors = c("#e06663", "#327eba"),
#                        guide = guide_colorbar(reverse = TRUE,
#                                               order = 1))

# transcription factor family annotate
plot_data <- subset(tf_id_annotation, TF_ID %in% tf_id)
plot_data <- plot_data[order(plot_data$Family), ]
plot_data$TF_ID <- factor(plot_data$TF_ID, levels = plot_data$TF_ID)
tf_annot_plot <- ggplot(data = plot_data,
                        aes(x = TF_ID, y = 1, fill = Family)) +
  geom_tile() +
  ggsci::scale_fill_simpsons(alpha = .6) +
  ggfun::theme_nothing()

fig <- insert_top(tf_annot_plot, enrich_plot, height = 5) |>
  insert_bottom(enrich_pathway_plot, height = 50)
fig





ggsave(fig, file = "Phyllostachys_heterocycla_example/result/fig.png",
       width = 15, height = 10)
ggsave(fig, file = "Phyllostachys_heterocycla_example/result/fig.pdf",
       width = 15, height = 10)
