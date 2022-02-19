#!/usr/bin/env Rscript 

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))

# Load common functions

suppressPackageStartupMessages(require(workflowscriptscommon))

# parse options

option_list = list(
  make_option(
    c("-i", "--input-object-files"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Comman separated list of RDS/Loom/SCE serialised objects to integrate. They should all be of the same format."
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
    c("--reference"),
    action = "store",
    default = NULL,
    type = 'character',
    metavar = "Reference",
    help = "A command separated vector specifying the object/s to be used as a reference during integration. If NULL (default), all pairwise anchors are found (no reference/s). If not NULL, the corresponding objects in object.list will be used as references. When using a set of specified references, anchors are first found between each query and each reference. The references are then integrated through pairwise integration. Each query is then mapped to the integrated reference."
  ),
  make_option(
    c("--assay-list"),
    action = "store",
    default = NULL,
    type = 'character',
    metavar = "Assays for each provided dataset",
    help = "Comma separated vector of assay names specifying which assay to use when constructing anchors. If NULL, the current default assay for each object is used."
  ),
  make_option(
    c("--anchor-features"),
    action = "store",
    default = NULL,
    type = 'character',
    metavar = "Anchor features",
    help = "Either a numeric value stating the number of anchor features or comma separated ."
  ),
  make_option(
    c("-n", "--normalization-method"),
    action = "store",
    default = NULL,
    metavar = "Normalization method",
    type = 'character',
    help = "Name of normalization method used: either LogNormalize or SCT"
  ),
  make_option(
    c("-s", "--do-scale"),
    action = "store_true",
    default = FALSE,
    metavar = "Scale",
    type = 'logical',
    help = "Whether or not to scale the features provided. Only set to FALSE if you have previously scaled the features you want to use for each object in the object.list. False by default."
  ),
  make_option(
    c("--sct-clip-range"),
    action = "store",
    default = NULL,
    type = 'character',
    metavar = "SCT Clipe range",
    help = "Numeric of length two (comma separated) specifying the min and max values the Pearson residual will be clipped to"
  ),  
  make_option(
    c("--reduction"),
    action = "store",
    default = NULL,
    metavar = "Dimensional reduction",
    type = 'double',
    help = "Dimensional reduction to perform when finding anchors, either cca (Canonical correlation analysis) or rpca (Reciprocal PCA)."
  ),
  make_option(
    c("--l2-norm"),
    action = "store_true",
    default = FALSE,
    type = 'logical',
    help = "Perform L2 normalization on the CCA cell embeddings after dimensional reduction."
  ),
  make_option(
    c("-d", "--dims"),
    action = "store",
    default = '1:30',
    metavar = "Dimensions to use from CCA",
    type = 'character',
    help = "Which dimensions to use from the CCA to specify the neighbor search space. Comma or conlon separated (for ranges)."
  ),
  make_option(
    c("--k-anchor"),
    action = "store",
    default = 5,
    metavar = "Number of neighbours for picking anchors",
    type = 'integer',
    help = "How many neighbors (k) to use when picking anchors."
  ),
  make_option(
    c("--k-filter"),
    action = "store",
    default = 200,
    metavar = "Number of neighbours for filtering anchors",
    type = 'integer',
    help = "How many neighbors (k) to use when filtering anchors."
  ),
  make_option(
    c("--k-score"),
    action = "store",
    default = 30,
    metavar = "Number of neighbours for scoring anchors",
    type = 'integer',
    help = "How many neighbors (k) to use when scoring anchors"
  ),
  make_option(
    c("--k-weight"),
    action = "store",
    default = 100,
    metavar = "Neighbors for weighting",
    type = 'integer',
    help = "Number of neighbors to consider when weighting"
  ),
  make_option(
    c("--max-features"),
    action = "store",
    default = 200,
    metavar = "Max features to specify neighborhood during filtering",
    type = 'integer',
    help = "The maximum number of features to use when specifying the neighborhood search space in the anchor filtering"
  ),
  make_option(
    c("--nn-method"),
    action = "store",
    default = "rann",
    metavar = "Max features to specify neighborhood during filtering",
    type = 'character',
    help = "Method for nearest neighbor finding. Either rann or annoy"
  ),
  make_option(
    c("--eps"),
    action = "store",
    default = 0,
    metavar = "Error boundary",
    type = 'integer',
    help = "Error bound on the neighbor finding algorithm (from RANN)"
  ),
  make_option(
    c("--features"),
    action = "store",
    default = NULL,
    metavar = "Features for weighting PCA",
    type = 'character',
    help = "List of features (comma separated) to use when computing the PCA to determine the weights. Only set if you want a different set from those used in the anchor finding process"
  ),
  make_option(
    c("--features-to-integrate"),
    action = "store",
    default = NULL,
    metavar = "Features for weighting PCA",
    type = 'character',
    help = "List of features (comman separated) to integrate. By default, will use the features used in anchor finding."
  ),
  make_option(
    c("--weight-reduction"),
    action = "store",
    default = NULL,
    metavar = "Dimension reduction to use for anchor weights",
    type = 'character',
    help = "Dimension reduction to use when calculating anchor weights. This can be either the name of a dimension reduction present in all objects to be integrated, vector of strings, specifying the name of a dimension reduction to use for each object to be integrated or NULL, in which case a new PCA will be calculated and used to calculate anchor weights"
  ),
  make_option(
    c("--sd-weight"),
    action = "store",
    default = 1,
    metavar = "Weighting bandwidth",
    type = 'float',
    help = "Controls the bandwidth of the Gaussian kernel for weighting"
  ),
  make_option(
    c("--sample-tree"),
    action = "store",
    default = NULL,
    metavar = "Sample tree",
    type = 'character',
    help = "Specify the order of integration. If NULL, will compute automatically"
  ),
  make_option(
    c("--preserve-order"),
    action = "store_true",
    default = FALSE,
    metavar = "Preserve order on integration",
    type = 'logical',
    help = "Do not reorder objects based on size for each pairwise integration."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'Seurat'.'"
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file', 'vars_to_regress'))

# Check parameter values
inputs<-strsplit(opt$input_object_files,split = ",")
if ( lenght(inputs) <= 1 ) {
  stop("At least 2 input objects need to be given for integration running.")
}
for ( input in inputs ) {
  if ( ! file.exists(input)){
    stop((paste('Input file', input, 'does not exist.')))
  }
}

weight_reduction<-NULL
if (! is.null(opt$weight_reduction)) {
  if ( grepl(pattern=',', x=opt$weight_reduction ) ) {
    weight_reduction<-strsplit(opt$weight_reduction, split=",")
    if (length(weight_reduction) != length(inputs) ) {
      stop(cat("--weight-reduction and --input-object-files need to be of the same lenght, currently", length(weight_reduction),"and", length(inputs)))
    }
  } else {
    weight_reduction<-opt$weight_reduction
  }
}

expand_sequence <- function(dims_processed, dim_entry) {
  dims_processed[!dims_processed %in% c(dim_entry)]->dims_processed
  from<-unlist(strsplit(dim_entry, split=":"))[1]
  to<-unlist(strsplit(dim_entry, split=":"))[2]
  dims_processed<-append(dims_processed, seq(from=from, to=to, by=1))
}

dims_processed<-NULL
# handle dims listed as "1,2,3" or "1:3" or a combination of both.
if ( !is.null(opt$dims) ) {
  if ( grepl(pattern = ",", x=opt$dims)) {
    dims_processed<-unlist(strsplit(opt$dims, split = ","))
    for(dim_entry in dims_processed) {
      if ( grepl(pattern = ":", x=dim_entry)) {
        dims_processed<-expand_sequence(dims_processed, dim_entry)
      }
    }
    dims_processed<-as.numeric(dims_processed)
  } else if ( grepl(pattern = ":", x=opt$dims) ) {
    dims_processed<-expand_sequence(c(), opt$dims)
  } else {
    dims_processed<-opt$dims
  }
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
if(opt$input_format == "loom" | opt$output_format == "loom") {
  suppressPackageStartupMessages(require(loomR))
} else if(opt$input_format == "singlecellexperiment" | opt$output_format == "singlecellexperiment") {
  suppressPackageStartupMessages(require(scater))
}

# Input from serialized R object
objects_list<-list()
for (input in inputs) {
  seurat_object <- read_seurat3_object(input_path = opt$input_object_file, format = opt$input_format)  
  append(objects_list, seurat_object)->objects_list
}


anchors <- FindIntegrationAnchors(object.list = objects_list, 
                                  sct.clip.range = sct_clip_range, l2.norm = opt$l2_norm, 
                                  k.anchor = opt$k_anchor, 
                                  k.filter = opt$k_filter, 
                                  k.score = opt$k_score, 
                                  reduction = opt$reduction, 
                                  max.features = opt$max_features, 
                                  nn.method = opt$nn_method, 
                                  eps = opt$eps, 
                                  verbose = FALSE,
                                  assay = opt$assay_list, 
                                  reference = opt$reference,  
                                  normalization.method = opt$normalization_method,
                                  scale = opt$do_scale,
                                  anchor.features = opt$anchor_features,
                                  dims = dims_processsed)
integrated<-IntegrateData(anchorset = anchors, 
              normalization.method = opt$normalization_method,
              k.weight = opt$k_weight, 
              sd.weight = opt$sd_weight,
              features = strsplit(opt$features, split=","),
              features.to.integrate = strsplit(opt$features.to.integrate, split = ",") , 
              weight.reduction = weight_reduction, 
              sample.tree = opt$sample_tree, 
              preserve.order = opt$preserve_order, 
              eps = opt$eps,
              verbose = FALSE,
              dims = dims_processed)

# Output to a serialized R object
write_seurat3_object(seurat_object = integrated, 
                     output_path = opt$output_object_file, 
                     format = opt$output_format)
