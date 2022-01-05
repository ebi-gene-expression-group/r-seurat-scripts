#!/usr/bin/env Rscript
# This script has been automatically generated through
#
# YAML2RScript.py -i ../r-seurat-scripts/cli-generated/manually_crafted_YAML/seurat-integration.yaml -o ../r-seurat-scripts/seurat-integration.R
#
# to change this file edit the input YAML and re-run the above command

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

option_list <- list(
    make_option(
        c("-i", "--input-object-files"),
        action = "store",
        metavar = "Input files",
        type = "character",
        help = "A list of Seurat objects between which to find anchors for downstream integration. Input files in either RDS-Seurat, Loom or SCE, comma separated."
    ),
    make_option(
        c("--input-format"),
        action = "store",
        default = "seurat",
        metavar = "Input format",
        type = "character",
        help = "Either loom, seurat, anndata or singlecellexperiment for the input format to read."
    ),
    make_option(
        c("--reference-object-files"),
        action = "store",
        default = NULL,
        metavar = "Input files for reference objects",
        type = "character",
        help = "A vector specifying the object/s to be used as a reference during integration. If NULL (default), all pairwise anchors are found (no reference/s). If not NULL, the corresponding objects in this list will be used as references. When using a set of specified references, anchors are first found between each query and each reference. The references are then integrated through pairwise integration. Each query is then mapped to the integrated reference."
    ),
    make_option(
        c("--reference-format"),
        action = "store",
        default = "seurat",
        metavar = "Reference input format",
        type = "character",
        help = "Either loom, seurat, anndata or singlecellexperiment for the input format to read."
    ),
    make_option(
        c("--assay-list"),
        action = "store",
        default = NULL,
        type = "character",
        help = "A vector of assay names specifying which assay to use when constructing anchors. If NULL, the current default assay for each object is used."
    ),
    make_option(
        c("--anchor-features"),
        action = "store",
        default = 2000,
        type = "integer",
        help = "A numeric value (this will call 'SelectIntegrationFeatures' to select the provided number of features to be used in anchor finding) or a vector of features to be used as input to the anchor finding process (comma separated)"
    ),
    make_option(
        c("-s", "--do-not-scale"),
        action = "store_false",
        default = TRUE,
        type = "logical",
        help = "Whether or not to scale the features provided. Only call if you have previously scaled the features you want to use for each object in the object.list"
    ),
    make_option(
        c("-n", "--normalization-method"),
        action = "store",
        default = "LogNormalize",
        type = "character",
        help = "Name of normalization method used: LogNormalize or SCT"
    ),
    make_option(
        c("--sct-clip-range"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Numeric of length two specifying the min and max values the Pearson residual will be clipped to"
    ),
    make_option(
        c("--reduction"),
        action = "store",
        default = "cca",
        type = "character",
        help = "Dimensional reduction to perform when finding anchors. Can be cca (Canonical correlation analysis) or rpca (Reciprocal PCA)"
    ),
    make_option(
        c("--do-not-l2-norm"),
        action = "store_false",
        default = TRUE,
        type = "logical",
        help = "Perform L2 normalization on the CCA cell embeddings after dimensional reduction"
    ),
    make_option(
        c("-d", "--dims"),
        action = "store",
        default = "1:30",
        type = "character",
        help = "Which dimensions to use from the CCA to specify the neighbor search space"
    ),
    make_option(
        c("--k-anchor"),
        action = "store",
        default = 5,
        type = "integer",
        help = "How many neighbors (k) to use when picking anchors"
    ),
    make_option(
        c("--k-filter"),
        action = "store",
        default = 200,
        type = "integer",
        help = "How many neighbors (k) to use when filtering anchors"
    ),
    make_option(
        c("--k-score"),
        action = "store",
        default = 30,
        type = "integer",
        help = "How many neighbors (k) to use when scoring anchors"
    ),
    make_option(
        c("--max-features"),
        action = "store",
        default = 200,
        type = "integer",
        help = "The maximum number of features to use when specifying the neighborhood search space in the anchor filtering"
    ),
    make_option(
        c("--nn-method"),
        action = "store",
        default = "annoy",
        type = "character",
        help = "Method for nearest neighbor finding. Options include: rann, annoy"
    ),
    make_option(
        c("--n-trees"),
        action = "store",
        default = 50,
        type = "integer",
        help = "More trees gives higher precision when using annoy approximate nearest neighbor search"
    ),
    make_option(
        c("--eps"),
        action = "store",
        default = 0,
        type = "integer",
        help = "Error bound on the neighbor finding algorithm (from RANN)"
    ),
    make_option(
        c("--verbose"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Print progress bars and output"
    ),
    make_option(
        c("--new-assay-name"),
        action = "store",
        default = "integrated",
        type = "character",
        help = "Name for the new assay containing the integrated data"
    ),
    make_option(
        c("--integrate-features-pca"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Vector of features to use when computing the PCA to determine the weights. Only set if you want a different set from those used in the anchor finding process"
    ),
    make_option(
        c("--features-to-integrate"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Vector of features to integrate. By default, will use the features used in anchor finding."
    ),
    make_option(
        c("--integrate-dims"),
        action = "store",
        default = "1:30",
        type = "character",
        help = "Number of dimensions to use in the anchor weighting procedure"
    ),
    make_option(
        c("--k-weight"),
        action = "store",
        default = 100,
        type = "integer",
        help = "Number of neighbors to consider when weighting anchors"
    ),
    make_option(
        c("--weight-reduction"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Dimension reduction to use when calculating anchor weights. This can be one of: A string, specifying the name of a dimension reduction present in all objects to be integrated; A vector of strings, specifying the name of a dimension reduction to use for each object to be integrated; A vector of DimReduc objects, specifying the object to use for each object in the integration; NULL, in which case a new PCA will be calculated and used to calculate anchor weights. Note that, if specified, the requested dimension reduction will only be used for calculating anchor weights in the first merge between reference and query, as the merged object will subsequently contain more cells than was in query, and weights will need to be calculated for all cells in the object."
    ),
    make_option(
        c("--sd-weight"),
        action = "store",
        default = 1,
        type = "integer",
        help = "Controls the bandwidth of the Gaussian kernel for weighting"
    ),
    make_option(
        c("--sample-tree"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Specify the order of integration. If NULL, will compute automatically."
    ),
    make_option(
        c("--preserve-order"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Do not reorder objects based on size for each pairwise integration."
    ),
    make_option(
        c("--integrate-eps"),
        action = "store",
        default = 0,
        type = "integer",
        help = "Error bound on the neighbor finding algorithm (from 'RANN')"
    ),
    make_option(
        c("-o", "--output-object-file"),
        action = "store",
        type = "character",
        help = "FILE IN"
    ),
    make_option(
        c("--output-format"),
        action = "store",
        default = "seurat",
        type = "character",
        help = "Output format"
    )
)

opt <- wsc_parse_args(option_list, 
                      mandatory = c("input_object_files", "output_object_file"))
                # Check parameter values
inputs<-strsplit(opt$input_object_files,split = ",")[[1]]
if ( length(inputs) <= 1 ) {
  stop("At least 2 input objects need to be given for integration running.")
}
for ( input in inputs ) {
  if ( ! file.exists(input)){
    stop((paste('Input file', input, 'does not exist.')))
  }
}

references<-NULL
if (!is.null(opt$reference_objects)) {
  reference_paths<-strsplit(opt$reference_objects, split = ",")
  for (ref in reference_paths) {
    if ( ! file.exists(ref)) {
      stop((paste('Reference file', ref, 'does not exist.')))
    }
  }
  references<-reference_paths
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
  opt$dims<-dims_processed
}

dims_processed<-NULL
if( ! is.null(opt$integrate_dims) ){
  if ( grepl(pattern = ",", x=opt$integrate_dims)) {
    dims_processed<-unlist(strsplit(opt$integrate_dims, split = ","))
    for(dim_entry in dims_processed) {
      if ( grepl(pattern = ":", x=dim_entry)) {
        dims_processed<-expand_sequence(dims_processed, dim_entry)
      }
    }
    dims_processed<-as.numeric(dims_processed)
  } else if ( grepl(pattern = ":", x=opt$integrate_dims) ) {
    dims_processed<-expand_sequence(c(), opt$integrate_dims)
  } else {
    dims_processed<-opt$dims
  }
  opt$integrate_dims<-dims_processed
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

if ( !is.null(opt$anchor_features) && is.character(opt$anchor_features)) {
  # this could be a number as well.
  opt$anchor_features <- strsplit(opt$anchor_features, split = ",")
}
if(! is.null(opt$features_to_integrate)) {
  opt$features_to_integrate <- strsplit(opt$features_to_integrate, split = ",")
}

options(future.globals.maxSize = 8000 * 1024^2)

load_seurat4_packages_for_format(formats = c(opt$input_format, opt$reference_format))

seurat_objects <- read_multiple_seurat4_objects(input_path_list = opt$input_object_files,
                    format = opt$input_format)

reference_objects <- read_multiple_seurat4_objects(input_path_list = opt$reference_object_files,
                    format = opt$reference_format)

anchors <- FindIntegrationAnchors(object.list = seurat_objects,
                    assay = opt$assay_list,
                    reference = reference_objects,
                    anchor.features = opt$anchor_features,
                    scale = opt$do_not_scale,
                    normalization.method = opt$normalization_method,
                    sct.clip.range = opt$sct_clip_range,
                    reduction = opt$reduction,
                    l2.norm = opt$do_not_l2_norm,
                    dims = opt$dims,
                    k.anchor = opt$k_anchor,
                    k.filter = opt$k_filter,
                    k.score = opt$k_score,
                    max.features = opt$max_features,
                    nn.method = opt$nn_method,
                    n.trees = opt$n_trees,
                    eps = opt$eps,
                    verbose = opt$verbose)

integrated <- IntegrateData(anchorset = anchors,
                    new.assay.name = opt$new_assay_name,
                    normalization.method = opt$normalization_method,
                    features = opt$integrate_features_pca,
                    features.to.integrate = opt$features_to_integrate,
                    dims = opt$integrate_dims,
                    k.weight = opt$k_weight,
                    weight.reduction = opt$weight_reduction,
                    sd.weight = opt$sd_weight,
                    sample.tree = opt$sample_tree,
                    preserve.order = opt$preserve_order,
                    eps = opt$integrate_eps,
                    verbose = opt$verbose)

write_seurat4_object(seurat_object = integrated,
                    output_path = opt$output_object_file,
                    format = opt$output_format)
