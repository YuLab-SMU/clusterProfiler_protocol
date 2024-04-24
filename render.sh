#!/bin/bash

Rscript -e 'rmarkdown::render("index.Rmd")'
Rscript -e 'rmarkdown::render("IBD_index.Rmd", output_file = "index", output_dir = "examples/IBD_2_subtypes")'
Rscript -e 'rmarkdown::render("Phyllostachys_heterocycla_index.Rmd", output_file = "index", output_dir = "examples/Phyllostachys_heterocycla")'
Rscript -e 'rmarkdown::render("single_cell_index.Rmd", output_file = "index", output_dir = "examples/single_cell")'
Rscript -e 'rmarkdown::render("annot_data_index.Rmd", output_file = "index", output_dir = "examples/Phyllostachys_heterocycla/annot_data")'
Rscript -e 'rmarkdown::render("cell_marker_db_index.Rmd", output_file = "index", output_dir = "examples/single_cell/cell_marker_db")'
Rscript -e 'rmarkdown::render("filtered_gene_bc_matrices.Rmd", output_file = "index", output_dir = "examples/single_cell/filtered_gene_bc_matrices")'
Rscript -e 'rmarkdown::render("hg19_index.Rmd", output_file = "index", output_dir = "examples/single_cell/filtered_gene_bc_matrices/hg19")'