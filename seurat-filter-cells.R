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
    c("-s", "--subset-names"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Parameters to subset on. Eg, the name of a gene, PC1, a column name in object@meta.data, etc. Any argument that can be retreived using FetchData."
  ),
  make_option(
    c("-l", "--low-thresholds"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Low cutoffs for the parameters (default is -Inf)."
  ),
  make_option(
    c("-j", "--high-thresholds"),
    action = "store",
    default = NA,
    type = 'character',
    help = "High cutoffs for the parameters (default is Inf)."
  ),
  make_option(
    c("-c", "--cells-use"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Comma-separated list of cell names to use as a subset. Alternatively, text file with one cell per line."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'Seurat'.'"
  )
)

opt <- rsw_parse_args(option_list, mandatory = c('input_object_file', 'subset_names', 'output_object_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

# Input from serialized R object

seurat_object <- readRDS(opt$input_object_file)

# Are the metadata variables valid for this object?

subset_names <- split_string(opt$subset_names)
check_metadata(seurat_object, subset_names)

# Parse numeric fields

lt <- parse_numeric(opt, 'low_thresholds', -Inf, length(subset_names))
ht <- parse_numeric(opt, 'high_thresholds', Inf, length(subset_names))

# Check the cells_use

cells_use <- opt$cells_use
if (! is.null(cells_use)){
  if (file.exists(cells_use)){
    cells_use <- read.table("test_cells.txt", stringsAsFactors = FALSE)$V1
  }else{
    cells_use <- split_string(cells_use)
  }
  
  check_cells(cells_use, seurat_object)
}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))

filtered_seurat_object <- FilterCells(seurat_object, subset.names = subset_names, low.thresholds = lt, high.thresholds = ht, cells.use = cells_use)

# Output to a serialized R object

saveRDS(filtered_seurat_object, file = opt$output_object_file)
