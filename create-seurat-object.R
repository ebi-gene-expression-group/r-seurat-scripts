#!/usr/bin/env Rscript 

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))

# Source common functions

ca <- commandArgs()
script_dir <- dirname(sub('--file=', '', ca[grep('--file', ca)]))
source(file.path(script_dir, 'r-seurat-workflow-accessory.R'))

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
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'Seurat'. File will have extension as per the object type.'"
  )
)

opt <- rsw_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('Directory', opt$input_object_file, 'does not exist')))
}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))

# Input from serialized R object

sc_matrix <- readRDS(opt$input_object_file)

# Create the Seurat object

seurat_object <- CreateSeuratObject(sc_matrix)

# Output to a serialized R object

saveRDS(seurat_object, file = opt$output_object_file)
