<!-- README.md is generated from README.Rmd. Please edit that file -->

# Using clusterProfiler to characterise Multi-Omics Data

If you use this work in published research, please cite:

Using clusterProfiler to characterise Multi-Omics Data

This repo contains source code and data to produce Figures of the above
paper.

The `IBD_2_subtypes_example`, `Phyllostachys_heterocyla_example` and
`single_cell_example` contain the data, scripts and results of the three
examples in the above article. Each sub directory contains `input_data`,
`result`, `script`.

  - input\_data: contains all the data sets that used to generate the
    figures.
  - result: contains the results.
  - script: contains the source code to produce the figures.

More details information can be found from
[here](https://yulab-smu.top/clusterProfiler_protocol/).

## Dependencies and locations

Here is the output of `sessionInfo()` of the system was compiled:

    ## R version 4.3.0 (2023-04-21)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /mnt/d/UbuntuApps/R/4.3.0/lib/R/lib/libRblas.so 
    ## LAPACK: /mnt/d/UbuntuApps/R/4.3.0/lib/R/lib/libRlapack.so;  LAPACK version 3.11.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Asia/Shanghai
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3               patchwork_1.2.0            
    ##  [3] ggsc_1.1.4                  ggrepel_0.9.3              
    ##  [5] CelliD_1.8.1                SingleCellExperiment_1.22.0
    ##  [7] SeuratObject_4.1.3          Seurat_4.3.0               
    ##  [9] ggfun_0.1.5                 DESeq2_1.40.1              
    ## [11] SummarizedExperiment_1.30.1 Biobase_2.60.0             
    ## [13] MatrixGenerics_1.12.0       matrixStats_0.63.0         
    ## [15] GenomicRanges_1.52.0        GenomeInfoDb_1.36.0        
    ## [17] IRanges_2.36.0              S4Vectors_0.40.2           
    ## [19] BiocGenerics_0.48.1         aplot_0.2.2                
    ## [21] dplyr_1.1.2                 enrichplot_1.22.0          
    ## [23] ggplot2_3.5.0               clusterProfiler_4.10.1     
    ## [25] MicrobiotaProcess_1.15.0    tictoc_1.2.1               
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] fs_1.6.2                  spatstat.sparse_3.0-1    
    ##   [3] bitops_1.0-7              HDO.db_0.99.1            
    ##   [5] httr_1.4.5                RColorBrewer_1.1-3       
    ##   [7] tools_4.3.0               sctransform_0.3.5        
    ##   [9] utf8_1.2.3                R6_2.5.1                 
    ##  [11] vegan_2.6-4               lazyeval_0.2.2           
    ##  [13] uwot_0.1.14               mgcv_1.8-42              
    ##  [15] permute_0.9-7             withr_2.5.0              
    ##  [17] sp_1.6-0                  progressr_0.13.0         
    ##  [19] cli_3.6.1                 spatstat.explore_3.1-0   
    ##  [21] scatterpie_0.2.2          sandwich_3.0-2           
    ##  [23] mvtnorm_1.1-3             spatstat.data_3.0-1      
    ##  [25] askpass_1.1               ggridges_0.5.4           
    ##  [27] pbapply_1.7-0             yulab.utils_0.1.4        
    ##  [29] gson_0.1.0                DOSE_3.26.1              
    ##  [31] scater_1.28.0             parallelly_1.35.0        
    ##  [33] RSQLite_2.3.1             generics_0.1.3           
    ##  [35] gridGraphics_0.5-1        ica_1.0-3                
    ##  [37] spatstat.random_3.1-5     GO.db_3.17.0             
    ##  [39] Matrix_1.5-4              ggbeeswarm_0.7.2         
    ##  [41] fansi_1.0.4               abind_1.4-5              
    ##  [43] lifecycle_1.0.3           multcomp_1.4-25          
    ##  [45] yaml_2.3.7                qvalue_2.32.0            
    ##  [47] SparseArray_1.2.4         Rtsne_0.16               
    ##  [49] grid_4.3.0                blob_1.2.4               
    ##  [51] promises_1.2.0.1          crayon_1.5.2             
    ##  [53] miniUI_0.1.1.1            lattice_0.21-8           
    ##  [55] beachmat_2.19.1           cowplot_1.1.1            
    ##  [57] KEGGREST_1.40.0           pillar_1.9.0             
    ##  [59] knitr_1.43                fgsea_1.26.0             
    ##  [61] future.apply_1.10.0       codetools_0.2-19         
    ##  [63] fastmatch_1.1-3           leiden_0.4.3             
    ##  [65] glue_1.6.2                RcppArmadillo_0.12.2.0.0 
    ##  [67] data.table_1.14.8         vctrs_0.6.3              
    ##  [69] png_0.1-8                 treeio_1.27.0            
    ##  [71] gtable_0.3.3              cachem_1.0.8             
    ##  [73] xfun_0.39                 S4Arrays_1.3.3           
    ##  [75] mime_0.12                 libcoin_1.0-9            
    ##  [77] tidygraph_1.2.3           survival_3.5-5           
    ##  [79] iterators_1.0.14          ellipsis_0.3.2           
    ##  [81] fitdistrplus_1.1-11       TH.data_1.1-2            
    ##  [83] ROCR_1.0-11               nlme_3.1-162             
    ##  [85] ggtree_3.9.1              bit64_4.0.5              
    ##  [87] RcppAnnoy_0.0.20          irlba_2.3.5.1            
    ##  [89] vipor_0.4.5               KernSmooth_2.23-22       
    ##  [91] colorspace_2.1-0          DBI_1.1.3                
    ##  [93] tidyselect_1.2.0          bit_4.0.5                
    ##  [95] compiler_4.3.0            BiocNeighbors_1.18.0     
    ##  [97] DelayedArray_0.29.4       plotly_4.10.1            
    ##  [99] shadowtext_0.1.2          scales_1.3.0             
    ## [101] lmtest_0.9-40             stringr_1.5.0            
    ## [103] digest_0.6.33             goftest_1.2-3            
    ## [105] spatstat.utils_3.0-3      rmarkdown_2.22           
    ## [107] XVector_0.40.0            htmltools_0.5.5          
    ## [109] pkgconfig_2.0.3           umap_0.2.10.0            
    ## [111] sparseMatrixStats_1.12.0  fastmap_1.1.1            
    ## [113] rlang_1.1.1               htmlwidgets_1.6.2        
    ## [115] shiny_1.7.4               DelayedMatrixStats_1.22.0
    ## [117] farver_2.1.1              zoo_1.8-12               
    ## [119] jsonlite_1.8.7            BiocParallel_1.34.2      
    ## [121] GOSemSim_2.27.2           BiocSingular_1.16.0      
    ## [123] RCurl_1.98-1.12           magrittr_2.0.3           
    ## [125] modeltools_0.2-23         scuttle_1.10.1           
    ## [127] GenomeInfoDbData_1.2.10   ggplotify_0.1.0          
    ## [129] munsell_0.5.0             Rcpp_1.0.10              
    ## [131] ape_5.7-1                 ggnewscale_0.4.9         
    ## [133] viridis_0.6.2             reticulate_1.28          
    ## [135] stringi_1.7.12            ggstar_1.0.4.001         
    ## [137] ggraph_2.1.0              zlibbioc_1.46.0          
    ## [139] MASS_7.3-59               plyr_1.8.8               
    ## [141] parallel_4.3.0            listenv_0.9.0            
    ## [143] deldir_1.0-6              Biostrings_2.68.1        
    ## [145] graphlayouts_1.0.0        splines_4.3.0            
    ## [147] tensor_1.5                locfit_1.5-9.7           
    ## [149] igraph_1.4.2              spatstat.geom_3.2-1      
    ## [151] ggtreeExtra_1.11.0        ggsignif_0.6.4           
    ## [153] ScaledMatrix_1.8.1        reshape2_1.4.4           
    ## [155] evaluate_0.21             RcppParallel_5.1.7       
    ## [157] foreach_1.5.2             tweenr_2.0.2             
    ## [159] httpuv_1.6.11             openssl_2.0.6            
    ## [161] RANN_2.6.1                tidyr_1.3.0              
    ## [163] purrr_1.0.1               polyclip_1.10-4          
    ## [165] future_1.32.0             scattermore_0.8          
    ## [167] ggforce_0.4.1             rsvd_1.0.5               
    ## [169] coin_1.4-2                xtable_1.8-4             
    ## [171] RSpectra_0.16-1           tidytree_0.4.5           
    ## [173] tidydr_0.0.5              later_1.3.1              
    ## [175] viridisLite_0.4.2         tibble_3.2.1             
    ## [177] beeswarm_0.4.0            memoise_2.0.1            
    ## [179] AnnotationDbi_1.62.1      cluster_2.1.4            
    ## [181] globals_0.16.2
