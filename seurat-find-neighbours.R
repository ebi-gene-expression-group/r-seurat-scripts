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
    c("-d", "--distance-matrix"),
    action = "store_true",
    default = FALSE,
    help = "Boolean value of whether the provided matrix is a distance matrix; note, for objects of class dist, this parameter will be set automatically."
  ),
  make_option(
    c("-k", "--k-param"),
    action = "store",
    default = 20,
    type = 'integer',
    help = "Defines k for the k-nearest neighbor algorithm"
  ),
  make_option(
    c("--no-compute-snn"),
    action = "store_false",
    default = TRUE,
    help = "Avoid computing the shared nearest neighbor graph. Default is to have it computed."
  ),
  make_option(
    c("--prune-snn"),
    action = "store",
    default = 1/15,
    type = 'double',
    help = "Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 --- prune everything)."
  ),
  make_option(
    c("--nn-method"),
    action = "store",
    default = 'rann',
    type = 'character',
    help = "Method for nearest neighbor finding. Options include: rann (default), annoy"
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'Seurat'.'"
  ),
  make_option(
    c("-a", "--annoy-metric"),
    action = "store",
    default = "euclidean",
    type = 'character',
    help = "Distance metric for annoy. Options include: euclidean (default), cosine, manhattan, and hamming"
  ),
  make_option(
    c("--graph-name"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Name of graph to use for the clustering algorithm."
  ),
  make_option(
    c("--nn-eps"),
    action = "store",
    default = 0.0,
    type = 'double',
    help = "Error bound when performing nearest neighbor seach using RANN; default of 0.0 implies exact nearest neighbor search"
  ),
  make_option(
    c("--verbose"),
    action = "store_true",
    default = FALSE,
    help = "Maximal number of iterations per random start"
  ),
  make_option(
    c("--force-recalc"),
    action = "store_true",
    default = FALSE,
    help = "Force recalculation of SNN"
  ),
  make_option(
    c("-f", "--features"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Comma-separated list of genes to use for building SNN. Alternatively, text file with one gene per line."
  ),
  make_option(
    c("--reduction"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Reduction to use as input for building the SNN"
  ),
  make_option(
    c("--dims"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Dimensions of reduction to use as input. A comma-separated list of the dimensions to use in construction of the SNN graph (e.g. To use the first 5 PCs, pass 1,2,3,4,5)."
  ),
  make_option(
    c("--assay"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Assay to use in construction of SNN"
  ),
  make_option(
    c("--do-plot"),
    action = "store_true",
    default = FALSE,
    help = "Plot SNN graph on tSNE coordinates"
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

dims <- opt$dims
if ( ! is.null(dims)){
  dims <- as.integer(wsc_parse_numeric(opt, 'dims'))
}

features <- NULL
if (! is.null(opt$features)){
  if (! file.exists(opt$features)){
    stop((paste('Supplied genes file', opt$features, 'does not exist')))
  }else{
    features <- readLines(opt$features)
  }
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

neighbours_object <- FindNeighbors(object = seurat_object,
                                  k.param = opt$k_param,
                                  compute.SNN = opt$no_compute_snn,
                                  prune.SNN = opt$prune_snn,
                                  nn.method = opt$nn_method,
                                  annoy.metric = opt$annoy_metric,
                                  nn.eps = opt$nn_eps,
                                  verbose = opt$verbose,
                                  force.recalc = opt$force_recalc,
                                  features = features,
                                  reduction = opt$reduction,
                                  dims = dims,
                                  assay = opt$assay,
                                  do.plot = opt$do_plot,
                                  graph.name = opt$graph_name
                                  )

# Summarise the clustering

# Some parameters aren't interesting for reporting purposes (e.g. file
# locations), so hide from the summary

nonreport_params <- c('input_object_file', 'output_object_file', 'help', 'output_text_file', 'tmp_file_location')
opt_table <- data.frame(value=unlist(opt), stringsAsFactors = FALSE)
opt_table <- opt_table[! rownames(opt_table) %in% nonreport_params, , drop = FALSE]

# Output to a serialized R object
write_seurat4_object(seurat_object = neighbours_object,
                     output_path = opt$output_object_file,
                     format = opt$output_format)
