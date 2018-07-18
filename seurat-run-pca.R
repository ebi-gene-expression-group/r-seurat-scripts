#!/usr/bin/env Rscript 

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))

# Source common functions

ca <- commandArgs()
script_dir <- dirname(sub('--file=', '', ca[grep('--file', ca)]))
source(file.path(script_dir, 'r-seurat-scripts-accessory.R'))

# parse options

option_list = list(
  make_option(
    c("-i", "--input-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which a serialized R matrix object may be found."
  ),
  make_option(
    c("-e", "--pc-genes"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "File to be used to derive a vector of gene names to scale/center. Default is all genes in object@data."
  ),
  make_option(
    c("-p", "--pcs-compute"),
    action = "store",
    default = 20,
    type = 'integer',
    help = "Total Number of PCs to compute and store (20 by default)."
  ),
  make_option(
    c("-m", "--use-imputed"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "Run PCA on imputed values (FALSE by default)."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'Seurat'.'"
  ),
  make_option(
    c("-b", "--output-embeddings-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store a csv-format embeddings table with PCs by cell."
  ),
  make_option(
    c("-l", "--output-loadings-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store a csv-format loadings table with PCs by gene."
  ),
  make_option(
    c("-s", "--output-stdev-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store PC stdev values (one per line)."
  )
)

opt <- rsw_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file', 'output_embeddings_file', 'output_loadings_file', 'output_stdev_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

if (! is.null(opt$pc_genes)){
  if (! file.exists(opt$pc_genes)){
    stop((paste('Supplied genes file', opt$genes_file, 'does not exist')))
  }else{
    pc_genes <- readLines(opt$pc_genes)
  }
}else{
  pc_genes <- NULL
}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))

# Input from serialized R object

seurat_object <- readRDS(opt$input_object_file)

pca_seurat_object <- RunPCA(seurat_object, pc.genes = pc_genes, pcs.compute = opt$pcs_compute, use.imputed = opt$use_imputed, do.print=FALSE)

# Output to text-format components

write.csv(pca_seurat_object@dr$pca@cell.embeddings, file = opt$output_embeddings_file)
write.csv(pca_seurat_object@dr$pca@gene.loadings, file = opt$output_loadings_file)
writeLines(con=opt$output_stdev_file, as.character(pca_seurat_object@dr$pca@sdev))

# Output to a serialized R object

saveRDS(pca_seurat_object, file = opt$output_object_file)

