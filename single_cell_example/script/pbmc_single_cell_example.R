library(Seurat)
library(CelliD)
library(clusterProfiler)
library(ggplot2)
library(patchwork)
# import cell marker database
cell_marker_db <- read.gmt(
    "./single_cell_example/input_data/cell_marker_db/c8.all.v2023.1.Hs.symbols.gmt"
)
# Load the PBMC dataset
pbmc.data <- Read10X(
    data.dir = "./single_cell_example/input_data/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
    min.cells = 3, min.features = 200)
# Preprocess
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc,
    subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# Normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize",
    scale.factor = 10000)
# Identification of highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunMCA(pbmc)
# Cluster the cells
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# cell type annotate using default method
seurat_cluster_id <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
    "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(seurat_cluster_id) <- levels(pbmc)
seurat_pbmc <- RenameIdents(pbmc, seurat_cluster_id)
seurat_pbmc_plot <- DimPlot(seurat_pbmc, label = TRUE, repel = TRUE) +
    theme(legend.position = "none")

# cell type annotate using clusterProfiler
predict_cell_type <- function(cell_type_enrich_result) {
    cell_type_enrich_result <- as.data.frame(cell_type_enrich_result)
    cell_type_enrich_result <- split(cell_type_enrich_result,
        f = cell_type_enrich_result$Cluster)
    sapply(cell_type_enrich_result, function(x) {
        x$ID[which.min(x$qvalue)]
    })
}
# Find cell group markers
marker_gene <- GetGroupGeneSet(X = pbmc, group.by = "seurat_clusters",
    n.features = 20)
cell_type_predict <- compareCluster(marker_gene, fun = "enricher",
    TERM2GENE = cell_marker_db) |> predict_cell_type()
clusterprofiler_pbmc <- RenameIdents(pbmc, cell_type_predict)
clusterprofiler_pbmc_plot <- DimPlot(
    clusterprofiler_pbmc, label = FALSE, repel = TRUE) +
    theme(legend.position = "none")
# compare two methods
seurat_pbmc_plot | clusterprofiler_pbmc_plot


