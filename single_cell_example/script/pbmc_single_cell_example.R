input_dir <- "single_cell_example/input_data"
library(tictoc)
tic("single cell example total")
tic("setup the environment and data objects")
library(Seurat)
library(CelliD)
library(clusterProfiler)
library(ggplot2)
library(ggrepel)
library(ggsc)
library(patchwork)

pbmc_counts <- Read10X(
  data.dir = file.path(input_dir, "filtered_gene_bc_matrices/hg19/")
)
pbmc <- CreateSeuratObject(counts = pbmc_counts, project = "pbmc3k",
                           min.cells = 3)
toc()

tic("Data preprocessing workflow")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize",
                      scale.factor = 10000)
pbmc <- ScaleData(pbmc)
toc()

tic("Dimensionality reduction")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",
                             nfeatures = 2000)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
toc()

tic("Cluster cells and identify markers of cell clusters")
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunMCA(pbmc)
cluster_markers <- GetGroupGeneSet(X = pbmc,
                                   group.by = "seurat_clusters",
                                   n.features = 20)
toc()

tic("cell type annotation")
cell_marker_db <- read.gmt(
  file.path(input_dir, "cell_marker_db/c8.all.v2023.1.Hs.symbols.gmt")
)
cell_type_enrich_result <- compareCluster(cluster_markers,
                                          fun = "enricher",
                                          TERM2GENE = cell_marker_db)

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

cols <- c('#B3D3AA', '#E88A71', '#DEAB76', '#CD574D',  
        '#BF38AE', '#176D84', '#7D83B7', '#4040C6', '#994B41')
clusterprofiler_pbmc_plot <- sc_dim(clusterprofiler_pbmc, geom = geom_point, size=.5) +
  sc_dim_geom_label(geom = ggrepel::geom_text_repel, 
                    color = "black", bg.color = "white") +
  scale_color_discrete(type=cols) +
  theme(legend.position = "none")
ggsave(clusterprofiler_pbmc_plot,
       file = "single_cell_example/result/figure4-cp.pdf",
       width = 12, height = 7)
toc()
toc()

# compare two methods
library(aplot)
# cell type annotate using default method
seurat_cluster_id <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B",
                       "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(seurat_cluster_id) <- levels(pbmc)
seurat_pbmc <- RenameIdents(pbmc, seurat_cluster_id)
seurat_pbmc_plot <- sc_dim(seurat_pbmc, geom = geom_point, size =.5) +
  sc_dim_geom_label(geom = ggrepel::geom_text_repel, color = "black",
                    bg.color = "white") +
  scale_color_discrete(type=cols) +
  theme(legend.position = "none")

d <- rbind(seurat_cluster_id, cell_type_predict)
rownames(d) <- c("Known cell type", "Predicted cell type")


library(gridExtra)
dd <- layer_data(clusterprofiler_pbmc_plot)
dd <- unique(dd[, c("colour", "group")])
tabfig <- tableGrob(
  d, theme = ttheme_default(base_size=10,
    colhead = list(bg_params = list(fill = dd$colour[order(dd$group)]))
  )
)

fig <- (seurat_pbmc_plot | clusterprofiler_pbmc_plot) / tabfig +
  plot_layout(heights = c(1, .2)) + plot_annotation(tag_levels = "A")


fig

ggsave(fig, file = "single_cell_example/result/figure4.pdf",
  width = 15.5, height = 7
)
ggsave(fig, file = "single_cell_example/result/figure4.png",
  width = 15.5, height = 7
)
