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
    c("-r", "--reduction-use"),
    action = "store",
    default = 'pca',
    type = 'character',
    help = 'Which dimensional reduction (e.g. PCA, ICA) to use for the tSNE. Default is PCA.'
  ),
  make_option(
    c("-c", "--cells-use"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "File to be used to derive a vector of which cells to analyze (default, all cells)."
  ),
  make_option(
    c("-d", "--dims-use"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "A comma-separated list of the which dimensions to use as input features."
  ),
  make_option(
    c("-e", "--genes-use"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "File to be used to derive a vector of gene names. If set, run the tSNE on this subset of genes (instead of running on a set of reduced dimensions). Not set (NULL) by default."
  ),
  make_option(
    c("-f", "--do-fast"),
    action = "store",
    default = TRUE,
    type = 'character',
    help = "If TRUE, uses the Barnes-hut implementation, which runs faster, but is less flexible. TRUE by default."
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
  )
)

opt <- rsw_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file', 'output_embeddings_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

# Check genes_use

if (! is.null(opt$genes_use) && opt$genes_use != 'NULL'){
  if (! file.exists(opt$genes_use)){
    stop((paste('Supplied genes file', opt$genes_use, 'does not exist')))
  }else{
    genes_use <- readLines(opt$genes_use)
  }
}else{
  genes_use <- NULL
}

# Read cells file (if present)

if (! is.null(opt$cells_use) && opt$cells_use != 'NULL'){
  if (! file.exists(opt$cells_use)){
    stop((paste('Supplied genes file', opt$cells_use, 'does not exist')))
  }else{
    cells_use <- readLines(opt$cells_use)
  }
}else{
  cells_use <- NULL
}

# Check dims-use

dims_use <- opt$dims_use
if ( ! is.null(dims_use)){
  dims_use <- as.numeric(unlist(strsplit(opt$dims_use, ',')))
}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))

# Input from serialized R object

seurat_object <- readRDS(opt$input_object_file)

tsne_seurat_object <- RunTSNE( seurat_object, reduction.use = opt$reduction_use, cells.use = cells_use, dims.use = dims_use, genes.use = NULL, do.fast = opt$do_fast  )

# Output to text-format components

write.csv(tsne_seurat_object@dr$tsne@cell.embeddings, file = opt$output_embeddings_file)

# Output to a serialized R object

saveRDS(tsne_seurat_object, file = opt$output_object_file)

