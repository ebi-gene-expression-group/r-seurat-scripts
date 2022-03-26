#!/usr/bin/env Rscript 

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))

# Load common functions

suppressPackageStartupMessages(require(workflowscriptscommon))

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
    c("--input-format"),
    action = "store",
    default = "seurat",
    type = 'character',
    help = "Either loom, seurat, anndata or singlecellexperiment for the input format to read."
  ),
  make_option(
    c("--output-format"),
    action = "store",
    default = "seurat",
    type = 'character',
    help = "Either loom, seurat, anndata or singlecellexperiment for the output format."
  ),
  make_option(
    c("-n", "--normalization-method"),
    action = "store",
    default = 'LogNormalize',
    type = 'character',
    help = "Method for normalization. Default is log-normalization (LogNormalize). Can be 'CLR' or 'RC' additionally."
  ),
  make_option(
    c("-s", "--scale-factor"),
    action = "store",
    default = 10000,
    type = 'integer',
    help = "Sets the scale factor for cell-level normalization."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'Seurat'.'"
  ),
  make_option(
    c("--margin"),
    action = "store",
    default = 1,
    type = 'integer',
    help = "If performing CLR normalization, normalize across features (1) or cells (2)."
  ),
  make_option(
    c("--block-size"),
    action = "store",
    default = NULL,
    type = 'integer',
    help = "How many cells should be run in each chunk, will try to split evenly across threads"
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))
if(opt$input_format == "loom" | opt$output_format == "loom") {
  suppressPackageStartupMessages(require(SeuratDisk))
} else if(opt$input_format == "singlecellexperiment" | opt$output_format == "singlecellexperiment") {
  suppressPackageStartupMessages(require(scater))
}

# Input from serialized R object

seurat_object <- read_seurat4_object(input_path = opt$input_object_file, format = opt$input_format)
normalised_seurat_object <- NormalizeData(seurat_object, 
                                          normalization.method = opt$normalization_method, 
                                          scale.factor = opt$scale_factor, 
                                          margin = opt$margin, 
                                          block.size = opt$block_size,
                                          verbose = FALSE)

# Output to a serialized R object
write_seurat4_object(seurat_object = normalised_seurat_object, 
                     output_path = opt$output_object_file, 
                     format = opt$output_format)
