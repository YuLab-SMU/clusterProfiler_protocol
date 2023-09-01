# import helper function
source("./Phyllostachys_heterocycla_example/script/Phyllostachys_heterocycla_helper_function.R")
set.seed(1234)
# prepare db
plant_tf_db <- read.table(
    "./Phyllostachys_heterocycla_example/input_data/annot_data/regulation_from_motif_CE_Phe.txt",
    sep = "\t", header = FALSE)
tf_id_annotation <- read.table(
    "./Phyllostachys_heterocycla_example/input_data/annot_data/Phe_TF_list.txt",
    header = TRUE, sep = "\t")
plant_tf_db <- plant_tf_db[, c(1, 3)]
colnames(plant_tf_db) <- c("TERM", "GENE")

# prepare counts
mat <- counts_prepare("./Phyllostachys_heterocycla_example/input_data/counts")

# data analysis
seven_day_result <- DESeq_analysis(mat,
    column = "group", test_level = "168h", ref_level = "0h")

one_day_result <- DESeq_analysis(mat,
    column = "group", test_level = "24h", ref_level = "0h")

two_hour_result <- DESeq_analysis(mat,
    column = "group", test_level = "2h", ref_level = "0h")

all_result <- list(
    prepare_gsea_gene_list(se = seven_day_result),
    prepare_gsea_gene_list(se = one_day_result),
    prepare_gsea_gene_list(se = two_hour_result))

names(all_result) <- c("168h_vs_0h", "24h_vs_0h", "2h_vs_0h")

# GSEA analysis
compare_enrich_result <- compareCluster(
    all_result, fun = "GSEA", minGSSize = 10, maxGSSize = 500,
    pvalueCutoff = .05, TERM2GENE = plant_tf_db)

enrich_plot <- dotplot(compare_enrich_result, includeAll = TRUE, showCategory = 25) +
    scale_color_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
    guide = guide_colorbar(reverse=TRUE, order=1)) +
    guides(size = guide_legend(override.aes=list(shape=1))) +
    theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
    panel.grid.major.x = element_blank(), strip.placement = 'outside') +
    coord_flip() + theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.text.y = element_text(size = 20),
        axis.title.y = element_blank())

# transcription factor function annotate
go_db <- read.table(
    file = "./Phyllostachys_heterocycla_example/input_data/annot_data/Phe_GO_annotation.txt",
    header = TRUE, sep = "\t", quote = "")
tf_id <- unique(enrich_plot$data$ID)
enrich_pathway_plot <- tf_go_annot(compare_enrich_result, go_db) |>
    # enrich_heatmap_plot()
    dotplot(by = 'count', color='qvalue', showCategory = 3) +
    theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 30, size=8)) +
    xlab(NULL)

# transcription factor family annotate
plot_data <- subset(tf_id_annotation, TF_ID %in% tf_id)
plot_data <- plot_data[order(plot_data$Family), ]
plot_data$TF_ID <- factor(plot_data$TF_ID, levels = plot_data$TF_ID)
tf_annot_plot <- ggplot(
    data = plot_data, aes(x = TF_ID, y = "type", fill = Family)) +
    geom_tile() + scale_fill_manual(values = rainbow(11, alpha = .4)) + 
    ggfun::theme_nothing()

fig <- insert_top(tf_annot_plot, enrich_plot, height = 10) |>
    insert_bottom(enrich_pathway_plot, height = 50)

