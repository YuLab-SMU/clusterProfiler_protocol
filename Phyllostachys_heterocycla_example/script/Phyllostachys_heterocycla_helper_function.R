library(tidyr)
library(DESeq2)
library(SummarizedExperiment)
library(clusterProfiler)
library(ggplot2)
library(aplot)

# helper function definition
counts_prepare <- function(count_files_path) {
    count_files <- list.files(count_files_path, full.names = TRUE)
    names(count_files) <- gsub(
        pattern = "\\.txt$",
        replacement = "", basename(count_files)
    )
    count_list <- lapply(count_files, function(count_file) {
        read.table(count_file, header = FALSE, row.names = 1) |>
            as.matrix()
    })
    mat <- do.call(cbind, count_list)
    colnames(mat) <- gsub(
        pattern = "counts(\\d+)\\-(\\d+)",
        replacement = "\\1h_rep\\2", names(count_list)
    )
    group_info <- data.frame(
        ID = colnames(mat),
        group = gsub(
            pattern = "_rep\\d+$",
            replacement = "", colnames(mat)
        )
    )
    gene <- gsub(
        pattern = "\\.EXON\\d+$",
        replacement = "", rownames(mat)
    )
    mat <- aggregate(mat, list(gene = gene), sum) |>
        subset(!grepl(pattern = "^__.+", x = gene))
    rownames(mat) <- mat$gene
    mat <- as.matrix(mat[, -1])
    count2se(counts = mat, group_info = group_info)
}

count2se <- function(
    counts, group_info, gene_info = NULL,
    gene_length = NULL) {
    if (is.null(gene_info)) {
        gene_info <- data.frame(gene_id = rownames(counts))
        rownames(gene_info) <- rownames(counts)
    }
    SummarizedExperiment(
        assays = list(
            counts = as.matrix(counts[
                rownames(gene_info),
                group_info$ID
            ])
        ),
        rowData = gene_info,
        colData = group_info
    )
}

DESeq_analysis <- function(se, column, ref_level = NULL,
    test_level = NULL) {
    dds <- DESeqDataSet(se,
        design = formula(paste(c("~", column), collapse = " ")))
    dds[[column]] <- relevel(dds[[column]], ref = ref_level)
    res <- DESeq(dds) |> results(tidy = TRUE,
        contrast = c(column, test_level, ref_level))
    rownames(res) <- res$row
    res <- res[rownames(se), ]
    rowData(se) <- cbind(rowData(se), res[, -1])
    se
}

prepare_gsea_gene_list <- function(se, keytype = NULL) {
    res <- rowData(se)
    if (is.null(keytype)) {
        genelist <- setNames(object = res$log2FoldChange,
            nm = rownames(res))
    } else {
        genelist <- setNames(object = res$log2FoldChange,
            nm = rowData(se)[[keytype]])
    }
    genelist <- genelist[!is.na(genelist)]
    genelist <- genelist[order(genelist, decreasing = TRUE)]
    genelist
}

gsea_analysis <- function(se, keytype = NULL,
    minGSSize = 10, maxGSSize = 500, pvalueCutoff = .05,
    TERM2GENE = NULL, TERM2NAME = NULL, gson = NULL, method = NULL, ...) {
    gene_list <- prepare_gsea_gene_list(se, keytype = keytype)
    if (!is.null(TERM2GENE)) {
        GSEA(geneList = gene_list, minGSSize = minGSSize,
            maxGSSize = maxGSSize, pvalueCutoff = pvalueCutoff,
            TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME)
    } else if (!is.null(gson)) {
        GSEA(geneList = gene_list, minGSSize = minGSSize,
            maxGSSize = maxGSSize, pvalueCutoff = pvalueCutoff,
            gson = gson)
    }
}

tf_go_annot <- function(enrich_result, go_db) {
    enrich_result_info <- as.data.frame(enrich_result)
    gene_list <- strsplit(enrich_result_info$core_enrichment, split = "/")
    names(gene_list) <- enrich_result_info$ID
    gene_list <- stack(gene_list)
    gene_list <- split(gene_list$values, gene_list$ind)
    compareCluster(gene_list, fun = "enricher",
        TERM2GENE = go_db[, c("GO_ID", "Gene_id")],
        TERM2NAME = go_db[, c("GO_ID", "GO_term")])
}

enrich_heatmap_plot <- function(enrich_result, label_format = 60) {
    default_labeller <- enrichplot:::default_labeller
    plot_data <- as.data.frame(enrich_result)
    plot_data$log_p <- -10 * log10(plot_data$qvalue)
    ggplot(data = plot_data, aes(x = Cluster, y = Description,
        fill = log_p)) + geom_raster() +
        scale_y_discrete(labels = default_labeller(label_format)) +
        scale_fill_gradient2(
            low = "#b3eebe", high = "#371ea3", mid = "#46bac2") +
        theme_bw() + theme_classic() + theme(axis.text.x = element_text(colour = "black",
        size = 6, vjust = 1, hjust = 1, angle = 30),
        axis.text.y = element_text(colour = "black",
        size = 6, hjust = 1), axis.title = element_blank())
}