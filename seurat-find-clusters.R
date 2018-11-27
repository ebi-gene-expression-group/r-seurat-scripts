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
    c("-e", "--genes-use"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "File to be used to derive a vector of gene names to use in construction of SNN graph if building directly based on expression data rather than a dimensionally reduced representation (i.e. PCs)."
  ),
  make_option(
    c("-u", "--reduction-type"),
    action = "store",
    default = 'pca',
    type = 'character',
    help = "Name of dimensional reduction technique to use in construction of SNN graph. (e.g. 'pca', 'ica')."
  ),
  make_option(
    c("-d", "--dims-use"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "A comma-separated list of the dimensions to use in construction of the SNN graph (e.g. To use the first 5 PCs, pass 1,2,3,4,5)."
  ),
  make_option(
    c("-k", "--k-param"),
    action = "store",
    default = 30,
    type = 'integer',
    help = "Defines k for the k-nearest neighbor algorithm."
  ),
  make_option(
    c("-j", "--prune-snn"),
    action = "store",
    default = 1/15,
    type = 'double',
    help = "Sets the cutoff for acceptable Jaccard distances when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Essentially sets the strigency of pruning (0 — no pruning, 1 — prune everything)."
  ),
  make_option(
    c("-r", "--resolution"),
    action = "store",
    default = 0.8,
    type = 'double',
    help = "Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."
  ),
  make_option(
    c("-a", "--algorithm"),
    action = "store",
    default = 1,
    type = 'integer',
    help = "Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm)."
  ),
  make_option(
    c("-m", "--tmp-file-location"),
    action = "store",
    default = 1,
    type = 'character',
    help = "Directory where intermediate files will be written. Specify the ABSOLUTE path."
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
    help = "File name in which to store text format set of clusters."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file', 'output_text_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

dims_use <- opt$dims_use
if ( ! is.null(dims_use)){
  dims_use <- wsc_parse_numeric(opt, 'dims_use')
}

if (! is.null(opt$genes_use)){
  if (! file.exists(opt$genes_use)){
    stop((paste('Supplied genes file', opt$genes_use, 'does not exist')))
  }else{
    genes_use <- readLines(opt$genes_use)
  }
}else{
  genes_use <- NULL
}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))

# Input from serialized R object

seurat_object <- readRDS(opt$input_object_file)

clustered_object <- FindClusters(seurat_object, genes.use = genes_use, reduction.type = opt$reduction_type, dims.use = dims_use, k.param = opt$k_param, prune.SNN = opt$prune_snn, print.output = FALSE, save.SNN = FALSE, resolution = opt$resolution, temp.file.location = opt$temp_file_location)

# Output to a serialized R object

saveRDS(clustered_object, file = opt$output_object_file)

# Output variable genes to a simple text file

write.csv(data.frame(cell=names(clustered_object@ident), cluster=clustered_object@ident), file = opt$output_text_file, row.names = FALSE)
