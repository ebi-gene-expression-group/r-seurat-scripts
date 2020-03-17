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
    help = "seurat object to use as the reference",
  ),
 
   make_option(
    c("--query-format"),
    action = "store",
    default = "seurat",
    type = 'character',
    help = "format of the query"
  ),
  make_option(
    c("-r", "--reference-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "seurat object to use as the reference"
  ),
  
  make_option(
    c("--reference-format"),
    action = "store",
    default = "seurat",
    type = 'character',
    help = "format of the reference"
  ),
  
  
  make_option(
    c("-o", "--output-file"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "path to the output"
  ),
  make_option(
    c("-n", "--normalization-method"),
    action = "store",
    default = 'SCT',
    type = 'character',
    help = "default SCT can also be LogNormalize"
  ),
  make_option(
    c("--I2-norm"),
    action = "store_false",
    default = FALSE,
    help = " execute a i2 normalization on the query"

  ),
  make_option(
    c("--reference-assay"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "assay to use from reference"
  ),
  make_option(
    c("--query-assay"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "assay to use from query"
  ),
  make_option(
    c("--reduction"),
    action = "store",
    default = "pcaproject",
    type = 'character',
    help = "Dimensional reduction to perform when finding anchors"
  ),
  make_option(
    c("--project-query"),
    action = "store_false",
    help = "Project the PCA from the query dataset onto the reference. Use only in rare cases"
  ),
  make_option(
    c("-f","--features"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Features to use for dimensional reductionFeatures to use for dimensional reduction"
  ),
  make_option(
    c("--npcs"),
    action = "store",
    default = 30,
    type = 'integer',
    help = "Number of PCs to compute on reference. If null, then use an existing PCA structure in the reference object"
  ),
  make_option(
    c("-d","--dims"),
    action = "store",
    default = 0:30,
    type = 'integer', #note sure
    help = "dimension to go throw"
  ),
  make_option(
    c("--k-anchor"),
    action = "store",
    default = 5,
    type = 'integer',
    help = "How many neighbors (k) to use when picking anchors"
  ),
  make_option(
    c("--k-filter"),
    
    action = "store",
    default = 200,
    type = 'integer',
    help = "How many neighbors (k) to use when filtering anchors"
  ),
  make_option(
    c("--k-score"),
    action = "store",
    default = 30,
    type = 'integer',
    help ="How many neighbors (k) to use when scoring anchors"
  ),
  make_option(
    c("-m","--max-features"),
    action = "store",
    default = 200,
    type = 'integer',
    help ="The maximum number of features to use when specifying the neighborhood search space in the anchor filtering"
  ),
  make_option(
    c("--nn-method"),
    action = "store",
    default = "rann",
    type = 'character',
    help ="Method for nearest neighbor finding. Options include: rann, annoy"
  ),
  make_option(
    c("--eps"),
    action = "store",
    default = 0,
    type = 'integer',
    help ="Error bound on the neighbor finding algorithm (from RANN)"
  ),
  make_option(
    c("--approx-pca"),
    action = "store_true",
    default = TRUE,
    help ="Use truncated singular value decomposition to approximate PCA"
  ),
  make_option(
    c("--verbose"),
    action = "store_true",
    default = TRUE,
    help ="Print progress bars and output"
  )
  
)

#minimum arguments to work
opt <- wsc_parse_args(option_list, mandatory = c('query_file', 'reference_file','output_file'))

# Check parameter values

if ( ! file.exists(opt$reference_file)){
  stop((paste('File', opt$reference_file, 'does not exist')))
}


#load seurat and packages needed to read input
suppressPackageStartupMessages(require(Seurat))


#if(opt$query_format == "loom" | opt$output_format == "loom") {
 # suppressPackageStartupMessages(require(loomR))
#} else if(opt$query_format == "singlecellexperiment" | opt$output_format == "singlecellexperiment") {
 # suppressPackageStartupMessages(require(scater))
#}

seurat_query <- read_seurat3_object(input_path = opt$reference_file, format = opt$reference_format)
seurat_reference <- read_seurat3_object(input_path = opt$query_file, format = opt$query_format)

#make the fonction work
anchor_object <- FindTransferAnchors(seurat_reference,
                                    seurat_query,
                                    normalization.method = opt$normalization.method,
                                    reference.assay = opt$reference.assay,
                                    query.assay = opt$query.assay,
                                    reduction = opt$reduction,
                                   # project.query = opt$project.query,
                                    features = opt$features,
                                    npcs = opt$npcs,
                                    l2.norm = opt$l2.norm,
                                    dims = opt$dims,
                                    k.anchor = opt$k.anchor,
                                    k.filter = opt$k.filter,
                                    k.score = opt$k.score,
                                    max.features = opt$max.features,
                                    nn.method = opt$nn.method,
                                    eps = opt$eps,
                                    approx.pca = opt$approx.pca,
                                    verbose = opt$verbose)

#directly save the anchorset
saveRDS(anchor_object, file = opt$output_file)


