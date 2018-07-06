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
    c("-m", "--mean-function"),
    action = "store",
    default = 'ExpMean',
    type = 'character',
    help = "Function to compute x-axis value (average expression). Default is to take the mean of the detected (i.e. non-zero) values."
  ),
  make_option(
    c("-d", "--dispersion-function"),
    action = "store",
    default = 'LogVMR',
    type = 'character',
    help = "Function to compute y-axis value (dispersion). Default is to take the standard deviation of all values."
  ),
  make_option(
    c("-l", "--x-low-cutoff"),
    action = "store",
    default = 0.1,
    type = 'integer',
    help = "Bottom cutoff on x-axis for identifying variable genes."
  ),
  make_option(
    c("-h", "--x-high-cutoff"),
    action = "store",
    default = 8,
    type = 'integer',
    help = "Top cutoff on x-axis for identifying variable genes."
  ),
  make_option(
    c("-l", "--y-low-cutoff"),
    action = "store",
    default = 1,
    type = 'integer',
    help = "Bottom cutoff on y-axis for identifying variable genes."
  ),
  make_option(
    c("-h", "--y-high-cutoff"),
    action = "store",
    default = Inf,
    type = 'integer',
    help = "Top cutoff on y-axis for identifying variable genes."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'Seurat'.'"
  ),
  make_option(
    c("-t", "--output-text-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store variable genes in plain text."
  )
)

opt <- rsw_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))

# Input from serialized R object

seurat_object <- readRDS(opt$input_object_file)
variable_genes_seurat_object <- FindVariableGenes(seurat_object, mean.function = get(opt$mean_function), dispersion.function = get(opt$dispersion_function), x.low.cutoff = opt$x_low_cutoff, x.high.cutoff = opt$x_high_cutoff, y.cutoff = opt$y_low_cutoff, y.high.cutoff = opt$y_high_cutoff, do.plot = FALSE, display.progress = FALSE)

# Output to a serialized R object

saveRDS(normalised_seurat_object, file = opt$output_object_file)

# Output variable genes to a simple text file

writeLines(con="~/Desktop/test.txt", variable_genes_seurat_object@var.genes)
