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
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))

# Input from serialized R object

seurat_object <- readRDS(opt$input_object_file)
normalised_seurat_object <- NormalizeData(seurat_object, 
                                          normalization.method = opt$normalization_method, 
                                          scale.factor = opt$scale_factor, 
                                          verbose = FALSE)

# Output to a serialized R object

saveRDS(normalised_seurat_object, file = opt$output_object_file)
