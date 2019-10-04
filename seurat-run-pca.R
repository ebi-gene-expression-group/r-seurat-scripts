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
    c("-e", "--pc-genes"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "File with gene names to scale/center. Default is all genes in object@data."
  ),
  make_option(
    c("-c", "--pc-cells"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "File with cell names to scale/center. Default is all cells in object@data."
  ),
  make_option(
    c("-p", "--pcs-compute"),
    action = "store",
    default = 20,
    type = 'integer',
    help = "Total Number of PCs to compute and store (20 by default)."
  ),
  make_option(
    c("-r", "--reverse-pca"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "Run PCA on reverse matrix (gene x cell; FALSE by default means cell x gene)."
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
  ),
  make_option(
    c("-l", "--output-loadings-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store a csv-format loadings table with PCs by gene."
  ),
  make_option(
    c("-s", "--output-stdev-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store PC stdev values (one per line)."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file', 'output_embeddings_file', 'output_loadings_file', 'output_stdev_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

pc_genes <- NULL
if (! is.null(opt$pc_genes)){
  if (! file.exists(opt$pc_genes)){
    stop((paste('Supplied genes file', opt$pc_genes, 'does not exist')))
  }else{
    pc_genes <- readLines(opt$pc_genes)
  }
}

pc_cells <- NULL
if (! is.null(opt$pc_cells)){
  if (! file.exists(opt$pc_cells)){
    stop((paste('Supplied cells file', opt$pc_cells, 'does not exist')))
  }else{
    pc_cells <- readLines(opt$pc_cells)
  }
}


# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))
if(opt$input_format == "loom" | opt$output_format == "loom") {
  suppressPackageStartupMessages(require(loomR))
} else if(opt$input_format == "singlecellexperiment" | opt$output_format == "singlecellexperiment") {
  suppressPackageStartupMessages(require(scater))
}

# Input from serialized R object

seurat_object <- read_seurat3_object(input_path = opt$input_object_file, format = opt$input_format)

features<-pc_genes
if(opt$reverse_pca) {
  features<-pc_cells
}
pca_seurat_object <- RunPCA(seurat_object, 
                            features = features, 
                            npcs = opt$pcs_compute, 
                            rev.pca = opt$reverse_pca,
                            verbose=FALSE)

# Output to text-format components
# Review question: Do we need to revert this for the reverse PCA case?
write.csv(pca_seurat_object[['pca']]@cell.embeddings, file = opt$output_embeddings_file)
write.csv(pca_seurat_object[['pca']]@feature.loadings, file = opt$output_loadings_file)
writeLines(con=opt$output_stdev_file, as.character(pca_seurat_object[['pca']]@stdev))

# Output to a serialized R object
write_seurat3_object(seurat_object = pca_seurat_object, 
                     output_path = opt$output_object_file, 
                     format = opt$output_format)

