#!/usr/bin/env Rscript 

# Given an R seurat object, fetch a set of random genes. This needs to be random
# but reproducible so we can use it in predictable tests

cl <- commandArgs(trailingOnly = TRUE)

input_object_file <- cl[1]
output_text_file <- cl[2]
ngenes <- as.numeric(cl[3])

if (! file.exists(input_object_file)){
  stop(paste('File', input_object_file, "does not exist"))
}

# Read Seurat object

suppressPackageStartupMessages(require(Seurat))
seurat_object <- readRDS(input_object_file)

# Set the seed to make random genes reproducible

set.seed(42)
random_genes <- sample(rownames(seurat_object@data), ngenes)

# Write to file

writeLines(con=output_text_file, random_genes)
