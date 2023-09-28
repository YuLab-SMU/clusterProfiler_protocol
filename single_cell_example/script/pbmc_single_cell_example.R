library(Seurat)
library(CelliD)
library(clusterProfiler)
library(ggplot2)
library(ggrepel)
library(ggsc)
library(patchwork)

input_dir <- "single_cell_example/input_data"
# import cell marker database
cell_marker_db <- read.gmt(
  file.path(input_dir, "cell_marker_db/c8.all.v2023.1.Hs.symbols.gmt")
)
# Load the PBMC dataset
pbmc_data <- Read10X(
  data.dir = file.path(input_dir, "filtered_gene_bc_matrices/hg19/")
)
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc_data, project = "pbmc3k",
  min.cells = 3
)
# Preprocess
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
)
# Normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Identification of highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",
  nfeatures = 2000
)
# Cluster the cells
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Find cluster biomarkers
pbmc <- RunMCA(pbmc)
cluster_markers <- GetGroupGeneSet(X = pbmc,
  group.by = "seurat_clusters",
  n.features = 20
)

# cell type annotate using default method
seurat_cluster_id <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B",
  "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet"
)
names(seurat_cluster_id) <- levels(pbmc)
seurat_pbmc <- RenameIdents(pbmc, seurat_cluster_id)
seurat_pbmc_plot <- DimPlot(seurat_pbmc, label = TRUE, repel = TRUE) +
  theme(legend.position = "none")



# predict cell type using clusterProfiler
cell_type_enrich_result <- compareCluster(cluster_markers,
  fun = "enricher",
  TERM2GENE = cell_marker_db
)

# cell type annotate using clusterProfiler
predict_cell_type <- function(enrich_result) {
  enrich_result <- as.data.frame(enrich_result)
  result <- split(
    enrich_result, enrich_result$Cluster
  ) |>
    sapply(function(x) {
      x$ID[which.min(x$qvalue)]
    })
  cell_type <- gsub("_", " ", result) |>
    yulab.utils::str_wrap(18)
  names(cell_type) <- names(result)
  return(cell_type)
}

cell_type_predict <- predict_cell_type(cell_type_enrich_result)
clusterprofiler_pbmc <- RenameIdents(pbmc, cell_type_predict)

clusterprofiler_pbmc_plot <- sc_dim(clusterprofiler_pbmc) +
  sc_dim_geom_label(geom = geom_text_repel, color = "black",
                    bg.color = "white") +
  theme(legend.position = "none")

# compare two methods
library(aplot)
fig <- seurat_pbmc_plot | clusterprofiler_pbmc_plot +
  plot_annotation(tag_levels = "A")


ggsave(fig, file = "single_cell_example/result/DimPlot.pdf",
  width = 12, height = 7
)
ggsave(fig, file = "single_cell_example/result/DimPlot.png",
  width = 12, height = 7
)

d <- rbind(seurat_cluster_id, cell_type_predict)
rownames(d) <- c("Known cell type", "Predicted cell type")


library(gridExtra)
dd <- layer_data(clusterprofiler_pbmc_plot)
dd <- unique(dd[, c("colour", "group")])
tabfig <- tableGrob(d, theme=ttheme_default(base_size=10, colhead=list(bg_params = list(fill=dd$colour[order(dd$group)]))))

#grid::grid.draw(tabfig)

fig2 <- (seurat_pbmc_plot | clusterprofiler_pbmc_plot) / tabfig +
  plot_layout(heights = c(1, .2)) + plot_annotation(tag_levels = "A")


#fig2 <- aplot::plot_list(seurat_pbmc_plot, clusterprofiler_pbmc_plot, ggplotify::as.ggplot(tabfig), 
#  design="AABB\nCCCC", tag_levels = 'A', heights=c(1, .2))

fig2

ggsave(fig2, file = "single_cell_example/result/DimPlot2.pdf",
  width = 15.5, height = 7
)
ggsave(fig2, file = "single_cell_example/result/DimPlot2.png",
  width = 15.5, height = 7
)
