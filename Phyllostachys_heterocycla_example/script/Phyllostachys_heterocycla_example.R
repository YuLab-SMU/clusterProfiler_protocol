# import helper function
source("./Phyllostachys_heterocycla_example/script/Phyllostachys_heterocycla_helper_function.R")

# prepare db
plant_tf_db <- read.table(
    "./Phyllostachys_heterocycla_example/input_data/annot_data/regulation_from_motif_CE_Phe.txt",
    sep = "\t", header = FALSE)
plant_tf_db <- plant_tf_db[, c(1, 3)]
colnames(plant_tf_db) <- c("TERM", "GENE")

tf_id_annotation <- read.table(
    "./Phyllostachys_heterocycla_example/input_data/annot_data/Phe_TF_list.txt",
    header = TRUE, sep = "\t")

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
    aes(shape = I(15)) +
    scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe")) +
    guides(size = guide_legend(override.aes=list(shape=0))) +
    theme_minimal() +
    theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
    panel.grid.major.x = element_blank()) +
    coord_flip() + theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.text.y = element_text(size = 15),
        axis.title.y = element_blank())

# transcription factor function annotate
go_db <- read.table(
    file = "./Phyllostachys_heterocycla_example/input_data/annot_data/Phe_GO_annotation.txt",
    header = TRUE, sep = "\t", quote = "")


# tf_genes <- geneInCategory(compare_enrich_result) |> lapply(ls2df) |> do.call('rbind', args = _) |> with(data = _, split(value, category))
# y = compareCluster(tf_genes, fun='enricher', TERM2GENE = go_db[, c(2,1)], TERM2NAME= go_db[, c(2, 3)])


tf_id <- unique(compare_enrich_result[,'ID'])
tf_genes <- split(plant_tf_db$GENE, plant_tf_db$TERM)[tf_id]
y = compareCluster(tf_genes, fun='enricher', TERM2GENE = go_db[, c(2,1)], TERM2NAME= go_db[, c(2, 3)])

enrich_pathway_plot <-  dotplot(y, by = 'count', color='qvalue', showCategory = 3, label_format=40) +
    guides(size = guide_legend(override.aes=list(shape=1))) +
    theme_minimal() +
    theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 30, size=8)) +
    xlab(NULL) + ggtitle(NULL) + 
    scale_color_gradientn(colors = c("#e06663", "#327eba"), guide = guide_colorbar(reverse=TRUE, order=1))

 
# transcription factor family annotate
plot_data <- subset(tf_id_annotation, TF_ID %in% tf_id)
plot_data <- plot_data[order(plot_data$Family), ]
plot_data$TF_ID <- factor(plot_data$TF_ID, levels = plot_data$TF_ID)
tf_annot_plot <- ggplot(
    data = plot_data, aes(x = TF_ID, y = 1, fill = Family)) +
    # geom_tile() + scale_fill_manual(values = rainbow(11, alpha = .4)) + 
    geom_tile() + 
    ggsci::scale_fill_simpsons(alpha=.6) +
    #scale_fill_manual(values = c("#63b2ee", "#76da91", "#f8cb7f", "#f89588", "#7cd6cf", "#9192ab", "#7898e1", "#efa666", "#eddd86", "#9987ce", "#63b2ee")) +
    theme_minimal() +
    ggfun::theme_nothing()

fig <- insert_top(tf_annot_plot, enrich_plot, height = 5) |>
    insert_bottom(enrich_pathway_plot, height = 50)


ggsave(fig, file = './Phyllostachys_heterocycla_example/result/fig.png', width=16, height=10)
ggsave(fig, file = './Phyllostachys_heterocycla_example/result/fig.pdf', width=16, height=10)
