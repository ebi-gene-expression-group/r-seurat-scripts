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
    c("-d", "--data-dir"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Directory containing the matrix.mtx, genes.tsv, and barcodes.tsv files provided by 10X. A vector or named vector can be given in order to load several data directories. If a named vector is given, the cell barcode names will be prefixed with the name."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R matrix object."
  )
)

opt <- rsw_parse_args(option_list, mandatory = c('data_dir', 'output_object_file'))

# Check parameter values

if ( ! dir.exists(opt$data_dir)){
  stop((paste('Directory', opt$data_dir, 'does not exist')))
}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))

# Read the data

sc_matrix <- Read10X(data.dir = opt$data_dir)

# Output to a serialized R object

saveRDS(sc_matrix, file = opt$output_object_file)
