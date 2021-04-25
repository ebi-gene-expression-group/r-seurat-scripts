#!/usr/bin/env Rscript 

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))

# Load common functions

suppressPackageStartupMessages(require(workflowscriptscommon))

# parse options

option_list = list(
  make_option(
    c("-f", "--data-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A tab-separated file containing expression data."
  ),
  make_option(
    c("-d", "--data-dir"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Directory containing the matrix.mtx, genes.tsv, and barcodes.tsv files matching 10X conventions (overrides --data-file)."
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
  ),
  make_option(
    c("--gene-column"),
    action = "store",
    default = 2,
    type = 'integer',
    metavar = 'Gene name column',
    help = "Specify which column of genes.tsv or features.tsv to use for gene names; default is 2."
  ),
  make_option(
    c("--not-unique-features"),
    action = "store_true",
    default = FALSE,
    type = 'logical',
    metavar = 'Do not make features unique',
    help = "Do not make feature names unique (default FALSE - make them unique)."
  ),
  make_option(
    c("--project"),
    action = "store",
    default = "SeuratProject",
    type = 'character',
    metavar = 'Sets the project name for the Seurat object.',
    help = "Sets the project name for the Seurat object."
  ),
  make_option(
    c("--names-field"),
    action = "store",
    default = NULL,
    type = 'integer',
    metavar = 'Index for field with cells name',
    help = "For the initial identity class for each cell, choose this field for the cell's name. E.g. If your cells are named as BARCODE_CLUSTER_CELLTYPE in the input matrix, set names.field to 3 to set the initial identities to CELLTYPE."
  ),
  make_option(
    c("--names-delim"),
    action = "store",
    default = NULL,
    type = 'character',
    metavar = 'Delimiter field within cells name',
    help = "For the initial identity class for each cell, choose this delimiter from the cell's column name. E.g. If your cells are named as BARCODE-CLUSTER-CELLTYPE, set this to '-' to separate the cell name into its component parts for picking the relevant field."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('output_object_file'))

if (is.na(opt$data_file) && is.na(opt$data_dir)){
    stop("One of --data-file or data-dir must be supplied")
}

# Check parameter values

if (! is.na(opt$data_dir)){
    if ( ! dir.exists(opt$data_dir)){
      stop((paste('Directory', opt$data_dir, 'does not exist')))
    }
}else{
    if ( ! file.exists(opt$data_file)){
      stop((paste('File', opt$data_file, 'does not exist')))
    }
}

cell_metadata<-NULL
if ( ! is.null(opt$metadata) ) {
  cell_metadata<-read.table(opt$metadata, 
                            header = TRUE, sep="\t",
                            row.names = 1, check.names = FALSE, 
                            stringsAsFactors = FALSE)
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

if ( ! is.na(opt$data_dir)){
    sc_matrix <- Read10X(data.dir = opt$data_dir, 
                         unique.features = !opt$not_unique_features,
                         gene.column = opt$gene_column)

    # Use the default show method to print feedback
    printSpMatrix2(sc_matrix, note.dropping.colnames = FALSE, maxp = 500)

}else{
    sc_matrix <- read.table(opt$data_file)
    print(paste(nrow(sc_matrix), 'x', ncol(sc_matrix), 'matrix of class', class(sc_matrix)))
}

seurat_object <- CreateSeuratObject(sc_matrix,
                                    min.cells = opt$min_cells, 
                                    min.features = opt$min_features, 
                                    meta.data = cell_metadata, 
                                    project = opt$project, 
                                    names.field = opt$names_field,
                                    names.delim = opt$names_delim
                                    )

if(opt$output_format == "loom") {
  suppressPackageStartupMessages(require(SeuratDisk))
} else if(opt$output_format == "singlecellexperiment") {
  suppressPackageStartupMessages(require(scater))
}

# Fix for loom:
# https://github.com/mojaveazure/loomR/issues/36
if (opt$output_format == "loom" ) {
  seurat_object <- FindVariableFeatures(seurat_object, verbose = FALSE)
}

cat(c(
  '# Object summary', 
  capture.output(print(seurat_object)), 
  '\n# Metadata sample', 
  capture.output(head(seurat_object@meta.data))
), 
sep = '\n')

# Output to a serialized R object
# Output to a serialized R object
write_seurat4_object(seurat_object = seurat_object, 
                     output_path = opt$output_object_file, 
                     format = opt$output_format)
