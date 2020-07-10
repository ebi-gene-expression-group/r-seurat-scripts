#!/usr/bin/env Rscript

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
    c("-o", "--output-file"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "File name in which to store serialized R matrix object."
  ),
  make_option(
    c("-n", "--normalization-method"),
    action = "store",
    default = 'SCT',
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
    c("-d","--dims"),
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
    default = TRUE,
    help ="Print progress bars and output."
  )
)

#minimum arguments to work
opt <- wsc_parse_args(option_list, mandatory = c('query_file', 'reference_file','output_file'))

# Check parameter values
if ( ! file.exists(opt$reference_file)){
  stop((paste('File', opt$reference_file, 'does not exist')))
}

#convert dims from "a:b" to a real vector
dims <- wsc_parse_numeric(opt, 'dims')
#load seurat and packages needed to read input
suppressPackageStartupMessages(require(Seurat))

#load loomR or scater if needed
if(opt$query_format == "loom") {
  suppressPackageStartupMessages(require(loomR))
} else if(opt$query_format == "singlecellexperiment") {
  suppressPackageStartupMessages(require(scater))
}

seurat_query <- read_seurat3_object(input_path = opt$reference_file, format = opt$reference_format)
seurat_reference <- read_seurat3_object(input_path = opt$query_file, format = opt$query_format)
#make the function work
anchor_object <- FindTransferAnchors(seurat_reference,
                                    seurat_query,
                                    normalization.method = opt$normalization_method,
                                    reference.assay = opt$reference_assay,
                                    query.assay = opt$query_assay,
                                    reduction = opt$reduction,
                                    project.query = opt$project_query,
                                    features = opt$features,
                                    npcs = opt$npcs,
                                    l2.norm = opt$l2_norm,
                                    dims,
                                    k.anchor = opt$k_anchor,
                                    k.filter = opt$k_filter,
                                    k.score = opt$k_score,
                                    max.features = opt$max_features,
                                    nn.method = opt$nn_method,
                                    eps = opt$eps,
                                    approx.pca = opt$approx_pca,
                                    verbose = opt$verbose)

#directly save the anchorset
saveRDS(anchor_object, file = opt$output_file)


