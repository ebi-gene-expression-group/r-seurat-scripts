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
    c("-r", "--reduction-use"),
    action = "store",
    default = 'pca',
    type = 'character',
    help = 'Which dimensional reduction (e.g. PCA, ICA) to use for the tSNE. Default is PCA.'
  ),
  make_option(
    c("--tsne-method"),
    action = "store",
    default = 'Rtsne',
    type = 'character',
    help = 'Select the method to use to compute the tSNE. Available methods are: Rtsne, Flt-SNE'
  ),
  make_option(
    c("--perplexity"),
    action = "store",
    default = NULL,
    type = 'integer',
    help = 'Perplexity value for tSNE, if none is set, the default from seurat is used.'
  ),
  make_option(
    c("-c", "--cells-use"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "File to be used to derive a vector of which cells to analyze (default, all cells)."
  ),
  make_option(
    c("--dim_embed"),
    action = "store",
    default = 2,
    type = 'integer',
    help = "The dimensional space of the resulting tSNE embedding (default is 2). For example, set to 3 for a 3d tSNE"
  ),
  make_option(
    c("-d", "--dims-use"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "A comma-separated list of the which dimensions to use as input features."
  ),
  make_option(
    c("-e", "--genes-use"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "File to be used to derive a vector of gene names. If set, run the tSNE on this subset of genes (instead of running on a set of reduced dimensions). Not set (NULL) by default."
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
    c("--random-seed"),
    action = "store",
    default = NULL,
    type = 'integer',
    help = "Seed of the random number generator"
  ),
  make_option(
    c("--add-iter"),
    action = "store_true",
    default = FALSE,
    metavar = "Add iterations",
    type = 'logical',
    help = "If an existing tSNE has already been computed, uses the current tSNE to seed the algorithm and then adds additional iterations on top of this"
  ),
  make_option(
    c("--reduction-key"),
    action = "store",
    default = 'tSNE_',
    metavar = 'Reductio key',
    type = 'character',
    help = 'dimensional reduction key, specifies the string before the number for the dimension names. tSNE_ by default'
  ),
  make_option(
    c("--reduction-name"),
    action = "store",
    default = "tsne",
    metavar = "Reduction name",
    type = "character",
    help = "dimensional reduction name, specifies the position in the object$dr list. tsne by default"
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file', 'output_embeddings_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

# Check genes_use
genes_use <- NULL
if (! is.null(opt$genes_use) && opt$genes_use != 'NULL'){
  if (! file.exists(opt$genes_use)){
    stop((paste('Supplied genes file', opt$genes_use, 'does not exist')))
  }else{
    genes_use <- readLines(opt$genes_use)
  }
}

# Read cells file (if present)
cells_use <- NULL
if (! is.null(opt$cells_use) && opt$cells_use != 'NULL'){
  if (! file.exists(opt$cells_use)){
    stop((paste('Supplied genes file', opt$cells_use, 'does not exist')))
  }else{
    cells_use <- readLines(opt$cells_use)
  }
}

# Check dims-use

dims_use <- opt$dims_use
if ( ! is.null(dims_use)){
  dims_use <- as.integer(wsc_parse_numeric(opt, 'dims_use'))
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

tsne_seurat_object <- RunTSNE( seurat_object, 
                               reduction = opt$reduction_use, 
                               tsne.method = opt$tsne_method,
                               dim.embed = opt$dim_embed,
                               cells = cells_use, 
                               dims = dims_use, 
                               seed.use = opt$random_seed, 
                               add.iter = opt$add_iter,
                               perplexity = opt$perplexity,
                               reduction.key = opt$reduction_key, 
                               reduction.name = opt$reduction_name, 
                               features = genes_use
                              )

# Output to text-format components

write.csv(tsne_seurat_object[['tsne']]@cell.embeddings, file = opt$output_embeddings_file)

# Output to a serialized R object

saveRDS(tsne_seurat_object, file = opt$output_object_file)

