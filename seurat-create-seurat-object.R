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
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'Seurat'. File will have extension as per the object type.'"
  ),
  make_option(
    c("--min_cells"),
    action = "store",
    default = NA,
    type = 'integer',
    help = "Include genes with detected expression in at least this many cells."
  ),
  make_option(
    c("--min_genes"),
    action = "store",
    default = NA,
    type = 'integer',
    help = 'Include cells where at least this many genes are detected.'
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('Directory', opt$input_object_file, 'does not exist')))
}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))

# Input from serialized R object

sc_matrix <- readRDS(opt$input_object_file)

# Create the Seurat object
if( !is.na(opt$min_cells) && !is.na(opt$min_genes)) {
  seurat_object <- CreateSeuratObject(sc_matrix, 
                                      min.cells = opt$min_cells,
                                      min.genes = opt$min_genes)
} else if( is.na(opt$min_cells) && !is.na(opt$min_genes)) {
  seurat_object <- CreateSeuratObject(sc_matrix, 
                                      min.genes = opt$min_genes,)
} else if( !is.na(opt$min_cells) && is.na(opt$min_genes)) {
  seurat_object <- CreateSeuratObject(sc_matrix, 
                                      min.cells = opt$min_cells)
} else {
  seurat_object <- CreateSeuratObject(sc_matrix )
}

# Print an object summary

cat(c(
  '# Object summary', 
  capture.output(print(seurat_object)), 
  '\n# Metadata sample', 
  capture.output(head(seurat_object@meta.data))
), 
sep = '\n')

# Output to a serialized R object

saveRDS(seurat_object, file = opt$output_object_file)
