#!/usr/bin/env Rscript 

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))

# Load common functions

suppressPackageStartupMessages(require(workflowscriptscommon))

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
  ),
  make_option(
    c("--output-format"),
    action = "store",
    default = "seurat",
    type = 'character',
    help = "Either loom, seurat, anndata or singlecellexperiment for the output format."
  ),
  make_option(
    c("--metadata"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Path to a file with metdata for the cells, first column should be cell identifiers as used in the cells 10x file."
  ),
  make_option(
    c("--min-cells"),
    action = "store",
    default = 0,
    type = 'integer',
    help = "Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff."
  ),
  make_option(
    c("--min-features"),
    action = "store",
    default = 0,
    type = 'integer',
    help = "Include cells where at least this many features are detected."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('data_dir', 'output_object_file'))

# Check parameter values

if ( ! dir.exists(opt$data_dir)){
  stop((paste('Directory', opt$data_dir, 'does not exist')))
}

cell_metadata<-NULL
if ( ! is.null(opt$metadata) ) {
  cell_metadata<-read.table(opt$metadata,header = TRUE, sep="\t", row.names = 1)
  # vvv below is to avoid https://github.com/satijalab/seurat/issues/2310
  for ( name in colnames(cell_metadata)) {
    cell_metadata[[name]]<-gsub("^$", "N/A", trimws(cell_metadata[[name]]))
    cell_metadata[[name]][is.na(cell_metadata[[name]])]<-"N/A"
  }
} 

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(Matrix))

# Read the data

sc_matrix <- Read10X(data.dir = opt$data_dir)

# Use the default show method to print feedback
printSpMatrix2(sc_matrix, note.dropping.colnames = FALSE, maxp = 500)

seurat_object <- CreateSeuratObject(sc_matrix,
                                    min.cells = opt$min_cells, 
                                    min.features = opt$min_features, 
                                    meta.data = cell_metadata
                                    )

cat(c(
  '# Object summary', 
  capture.output(print(seurat_object)), 
  '\n# Metadata sample', 
  capture.output(head(seurat_object@meta.data))
), 
sep = '\n')

# Fix for loom:
# https://github.com/mojaveazure/loomR/issues/36
if (opt$output_format == "loom" ) {
  seurat_object <- FindVariableFeatures(seurat_object, verbose = FALSE)
}

# Output to a serialized R object
# Output to a serialized R object
write_seurat3_object(seurat_object = seurat_object, 
                     output_path = opt$output_object_file, 
                     format = opt$output_format)
