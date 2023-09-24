#!/usr/bin/env Rscript

# This deals with the functionality listed in https://satijalab.org/seurat/articles/integration_mapping.html
# Cell type classification using an integrated reference.

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

#set parsing options

option_list = list(
   make_option(
    c("-i", "--query-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File containing seurat object to use as the query.",
   ),
   make_option(
    c("--query-format"),
    action = "store",
    default = "seurat",
    type = 'character',
    help = "Either loom, seurat, anndata or singlecellexperiment for the output format."
  ),
  make_option(
    c("-r", "--reference-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File containing seurat object to use as the reference."
  ),
  make_option(
    c("--reference-format"),
    action = "store",
    default = "seurat",
    type = 'character',
    help = "Either loom, seurat, anndata or singlecellexperiment for the output format."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'Seurat'.'"
  ),
  make_option(
    c("--output-format"),
    action = "store",
    default = "seurat",
    type = 'character',
    help = "Either loom, seurat, anndata or singlecellexperiment for the output format."
  ),
  make_option(
    c("--output-anchorset-file"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "File name in which to store serialized R object for the anchorset."
  ),
  make_option(
    c("--output-anchorset-format"),
    action = "store",
    default = "seurat",
    type = 'character',
    help = "Either loom, seurat, anndata or singlecellexperiment for the output format."
  ),
  make_option(
    c("-n", "--normalization-method"),
    action = "store",
    default = 'LogNormalize',
    type = 'character',
    help = "Name of normalization method used: LogNormalize or SCT."
  ),
  make_option(
    c("--reference-assay"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Assay to use from reference."
  ),
  make_option(
    c("--query-assay"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Assay to use from query."
  ),
  make_option(
    c("--reduction"),
    action = "store",
    default = "pcaproject",
    type = 'character',
    help = "Dimensional reduction to perform when finding anchors."
  ),
  make_option(
    c("--reference_reduction"),
    action = "store",
    default = "pca",
    type = 'character',
    help = "Name of dimensional reduction to use from the reference if running the pcaproject workflow. Optionally enables reuse of precomputed reference dimensional reduction. If NULL (default), use a PCA computed on the reference object."
  ),
  make_option(
    c("--project-query"),
    action = "store_false",
    default = FALSE,
    help = "Project the PCA from the query dataset onto the reference. Use only in rare cases."
  ),
  make_option(
    c("--l2-norm"),
    action = "store_true",
    default = TRUE,
    help = "Execute a l2 normalization on the query."
  ),
  make_option(
    c("-f","--features"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Features to use for dimensional reductionFeatures to use for dimensional reduction."
  ),
  make_option(
    c("--npcs"),
    action = "store",
    default = 30,
    type = 'integer',
    help = "Number of PCs to compute on reference. If null, then use an existing PCA structure in the reference object."
  ),
  make_option(
    c("--dims"),
    action = "store",
    default = "1:30",
    type = 'character',
    help = "Which dimensions to use from the reduction to specify the neighbor search space."
  ),
  make_option(
    c("--k-anchor"),
    action = "store",
    default = 5,
    type = 'integer',
    help = "How many neighbors (k) to use when picking anchors."
  ),
  make_option(
    c("--k-filter"),
    action = "store",
    default = 200,
    type = 'integer',
    help = "How many neighbors (k) to use when filtering anchors."
  ),
  make_option(
    c("--k-score"),
    action = "store",
    default = 30,
    type = 'integer',
    help ="How many neighbors (k) to use when scoring anchors."
  ),
  make_option(
    c("-m","--max-features"),
    action = "store",
    default = 200,
    type = 'integer',
    help ="The maximum number of features to use when specifying the neighborhood search space in the anchor filtering."
  ),
  make_option(
    c("--nn-method"),
    action = "store",
    default = "rann",
    type = 'character',
    help ="Method for nearest neighbor finding. Options include: rann, annoy."
  ),
  make_option(
    c("--eps"),
    action = "store",
    default = 0,
    type = 'integer',
    help ="Error bound on the neighbor finding algorithm (from RANN)."
  ),
  make_option(
    c("--approx-pca"),
    action = "store_true",
    default = TRUE,
    help ="Use truncated singular value decomposition to approximate PCA."
  ),
  make_option(
    c("--verbose"),
    action = "store_true",
    default = FALSE,
    help ="Print progress bars and output."
  ),
  make_option(
    c("--transfer-refdata"),
    action = "store",
    type = 'character',
    help ="Data to transfer. The name of the metadata field or assay from the reference object provided. This requires the reference parameter to be specified. "
  ),
  make_option(
    c("--transfer-weight-reduction"),
    action = "store",
    default = 'pcaproject',
    type = 'character',
    help ="Dimensional reduction to use for the weighting anchors. Options are: pcaproject - Use the projected PCA used for anchor building; pca - Use an internal PCA on the query only; cca - Use the CCA used for anchor building; custom DimReduc - User provided DimReduc paths as RDS, computed on the query cells"
  ),
  make_option(
    c("--transfer-l2-norm"),
    action = "store_true",
    default = FALSE,
    help ="Perform L2 normalization on the cell embeddings after dimensional reduction"
  ),
  make_option(
    c("--transfer-dims"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Set of dimensions to use in the anchor weighting procedure."
  ),
  make_option(
    c("--transfer-k-weight"),
    action = "store",
    default = 50,
    type = 'integer',
    help = "Number of neighbors to consider when weighting anchors"
  ),
  make_option(
    c("--transfer-sd-weight"),
    action = "store",
    default = 1,
    type = 'integer',
    help = "Controls the bandwidth of the Gaussian kernel for weighting"
  ),
  make_option(
    c("--transfer-eps"),
    action = "store",
    default = 0,
    type = 'integer',
    help ="Error bound on the neighbor finding algorithm (from RANN), for transfer"
  ),
  make_option(
    c("--transfer-prediction-assay"),
    action = "store_true",
    default = FALSE,
    help = "Return an Assay object with the prediction scores for each class stored in the data slot."
  ),
  make_option(
    c("--transfer-n-trees"),
    action = "store",
    default = 50,
    type = 'integer',
    help = "More trees gives higher precision when using annoy approximate nearest neighbor search."
  ),
  make_option(
    c("--transfer-slot"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Slot to store the imputed data. Must be either 'data' (default) or 'counts'"
  ),
  make_option(
    c("--transfer-no-store-weights"),
    action = "store_false",
    default = TRUE,
    help ="Don't store the weights matrix used for predictions in the returned query object."
  ),
  make_option(
    c("--metadata-col"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Col name to store metadata in."
  )
)

#minimum arguments to work
opt <- wsc_parse_args(option_list, mandatory = c('query_file', 'reference_file','output_object_file', 'transfer_refdata', 'transfer_weight_reduction'))

# Check parameter values
if ( ! file.exists(opt$reference_file)){
  stop((paste('File', opt$reference_file, 'does not exist')))
}

#convert dims from "a:b" to a real vector
dims <- wsc_parse_numeric(opt, 'dims')
#load seurat and packages needed to read input
suppressPackageStartupMessages(require(Seurat))

#load loomR or scater if needed
if (opt$query_format == "loom" | opt$reference_format == "loom" | opt$output_format == "loom") {
  suppressPackageStartupMessages(require(SeuratDisk))
}
if (opt$query_format == "singlecellexperiment" | opt$reference_format == "singlecellexperiment" | opt$output_format == "singlecellexperiment") {
  suppressPackageStartupMessages(require(scater))
}

seurat_query <- read_seurat4_object(input_path = opt$query_file, format = opt$query_format)
seurat_reference <- read_seurat4_object(input_path = opt$reference_file, format = opt$reference_format)

if(!(opt$reference_reduction %in% names(seurat_reference@reductions))) {
  print(paste0("Calculating reduction ", opt$reference_reduction," as it is not present in the reference object."))
  if(opt$reference_reduction == "pca") {
    seurat_reference <- ScaleData(seurat_reference, verbose = FALSE)
    seurat_reference <- RunPCA(seurat_reference, npcs=30, verbose = FALSE)
  } else {
    stop((paste('Reduction default calculation on reference not implemented for ', opt$reference_reduction, ', please compute previously on reference object.')))
  }
}
#make the function work
anchor_object <- FindTransferAnchors(seurat_reference,
                                    seurat_query,
                                    normalization.method = opt$normalization_method,
                                    reference.assay = opt$reference_assay,
                                    query.assay = opt$query_assay,
                                    reduction = opt$reduction,
                                    reference.reduction = opt$reference_reduction,
                                    project.query = opt$project_query,
                                    features = opt$features,
                                    npcs = opt$npcs,
                                    l2.norm = opt$l2_norm,
                                    dims = dims,
                                    k.anchor = opt$k_anchor,
                                    k.filter = opt$k_filter,
                                    k.score = opt$k_score,
                                    max.features = opt$max_features,
                                    nn.method = opt$nn_method,
                                    eps = opt$eps,
                                    approx.pca = opt$approx_pca,
                                    verbose = opt$verbose)

#directly save the anchorset
if (!is.null(opt$output_anchorset_file)) {
  print(paste0("Output anchorset: ", opt$output_anchorset_file))
  write_seurat4_object(seurat_object = anchor_object,
                       output_path = opt$output_anchorset_file,
                       format = opt$output_anchorset_format
                       )
}

predictions<-TransferData(
  anchorset = anchor_object,
  refdata = opt$transfer_refdata,
  reference = seurat_reference,
  query = seurat_query,
  weight.reduction = opt$transfer_weight_reduction,
  l2.norm = opt$transfer_l2_norm,
  dims = opt$transfer_dims,
  k.weight = opt$transfer_k_weight,
  sd.weight = opt$transfer_sd_weight,
  eps = opt$transfer_eps,
  n.trees = opt$transfer_n_trees,
  verbose = opt$verbose,
  slot = opt$transfer_slot,
  prediction.assay = opt$transfer_prediction_assay,
  store.weights = opt$transfer_no_store_weights
)

seurat_query <- AddMetaData(seurat_query, metadata = predictions@meta.data, col.name = opt$metadata_col)

# Output to a serialized R object
write_seurat4_object(seurat_object = seurat_query,
                     output_path = opt$output_object_file,
                     format = opt$output_format)

