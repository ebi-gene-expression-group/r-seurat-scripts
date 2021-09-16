#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(SeuratDisk))
suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(workflowscriptscommon))

option_list <- list(
    make_option(
                    c('-i','--input-path'),
                    action='store',
                    metavar='Input file',
                    type='character',
                    help='Input file with Seurat object in either RDS-Seurat, Loom or SCE')
,
    make_option(
                    c('--input-format'),
                    action='store',
                    default='seurat',
                    metavar='Input format',
                    type='character',
                    help='Either loom, seurat, anndata or singlecellexperiment for the input format to read.')
,
    make_option(
                    c('--dims'),
                    action='store',
                    default=NULL,
                    type='character',
                    help='Which dimensions to use as input features, used only if list("features")  is NULL')
,
    make_option(
                    c('--reduction'),
                    action='store',
                    default='pca',
                    type='character',
                    help='Which dimensional reduction (PCA or ICA) to use for the UMAP input. Default is PCA')
,
    make_option(
                    c('--features'),
                    action='store',
                    default=NULL,
                    type='character',
                    help='Comma-separated list of features to use. If set, run UMAP on this subset of features (instead of running on a set of reduced dimensions). Not set (NULL) by default;  dims must be NULL to run on features')
,
    make_option(
                    c('--graph'),
                    action='store',
                    default=NULL,
                    type='character',
                    help='Name of graph on which to run UMAP')
,
    make_option(
                    c('--assay'),
                    action='store',
                    default=NULL,
                    type='character',
                    help='Assay to pull data for when using  list("features") , or assay used to construct Graph if running UMAP on a Graph')
,
    make_option(
                    c('--nn.name'),
                    action='store',
                    default=NULL,
                    type='character',
                    help='Name of knn output on which to run UMAP')
,
    make_option(
                    c('--slot'),
                    action='store',
                    default='data',
                    type='character',
                    help='The slot used to pull data for when using  list("features") . data slot is by default.')
,
    make_option(
                    c('--umap.method'),
                    action='store',
                    default='uwot',
                    type='character',
                    help='UMAP implementation to run. Can be list uwot, uwot-learn, umap-learn (rquires python umap-learn package).')
,
    make_option(
                    c('--reduction.model'),
                    action='store',
                    default=NULL,
                    type='character',
                    help='list("DimReduc")  object that contains the umap model')
,
    make_option(
                    c('--return.model'),
                    action='store_true',
                    default=FALSE,
                    type='logical',
                    help='whether UMAP will return the uwot model')
,
    make_option(
                    c('--n.neighbors'),
                    action='store',
                    default=30,
                    type='integer',
                    help='This determines the number of neighboring points used in local approximations of manifold structure. Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50.')
,
    make_option(
                    c('--n.components'),
                    action='store',
                    default=2,
                    type='integer',
                    help='The dimension of the space to embed into.')
,
    make_option(
                    c('--metric'),
                    action='store',
                    default='cosine',
                    type='character',
                    help='metric: This determines the choice of metric used to measure distance in the input space. A wide variety of metrics are already coded, and a user defined function can be passed as long as it has been JITd by numba.')
,
    make_option(
                    c('--n.epochs'),
                    action='store',
                    default=NULL,
                    type='integer',
                    help='The number of training epochs to be used in optimizing the low dimensional embedding. Larger values result in more accurate embeddings. If NULL is specified, a value will be selected based on the size of the input dataset (200 for large datasets, 500 for small).')
,
    make_option(
                    c('--learning.rate'),
                    action='store',
                    default=1,
                    type='integer',
                    help='The initial learning rate for the embedding optimization.')
,
    make_option(
                    c('--min.dist'),
                    action='store',
                    default=0,
                    type='integer',
                    help='This controls how tightly the embedding is allowed compress points together. Larger values ensure embedded points are moreevenly distributed, while smaller values allow the algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5.')
,
    make_option(
                    c('--spread'),
                    action='store',
                    default=1,
                    type='integer',
                    help='The effective scale of embedded points. In combination with min.dist this determines how clustered/clumped the embedded points are.')
,
    make_option(
                    c('--set.op.mix.ratio'),
                    action='store',
                    default=1.0,
                    type='double',
                    help='Interpolate between (fuzzy) union and intersection as the set operation used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial sets. Both fuzzy set operations use the product t-norm. The value of this parameter should be between 0.0 and 1.0; a value of 1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy intersection.')
,
    make_option(
                    c('--local.connectivity'),
                    action='store',
                    default=1,
                    type='integer',
                    help='The local connectivity required - i.e. the number of nearest neighbors that should be assumed to be connected at a local level. The higher this value the more connected the manifold becomes locally. In practice this should be not more than the local intrinsic dimension of the manifold.')
,
    make_option(
                    c('--repulsion.strength'),
                    action='store',
                    default=1,
                    type='integer',
                    help='Weighting applied to negative samples in low dimensional embedding optimization. Values higher than one will result in greater weight being given to negative samples.')
,
    make_option(
                    c('--negative.sample.rate'),
                    action='store',
                    default=5,
                    type='integer',
                    help='The number of negative samples to select per positive sample in the optimization process. Increasing this value will result in greater repulsive force being applied, greater optimization cost, but slightly more accuracy.')
,
    make_option(
                    c('--a'),
                    action='store',
                    default=NULL,
                    type='character',
                    help='More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread. Parameter of differentiable approximation of right adjoint functor.')
,
    make_option(
                    c('--b'),
                    action='store',
                    default=NULL,
                    type='character',
                    help='More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread. Parameter of differentiable approximation of right adjoint functor.')
,
    make_option(
                    c('--uwot.sgd'),
                    action='store_true',
                    default=FALSE,
                    type='logical',
                    help='Set list("uwot::umap(fast_sgd = TRUE)") ; see  list(list("umap"))  for more details')
,
    make_option(
                    c('--seed.use'),
                    action='store',
                    default=42,
                    type='integer',
                    help='Set a random seed. By default, sets the seed to 42. Setting NULL will not set a seed')
,
    make_option(
                    c('--metric.kwds'),
                    action='store',
                    default=NULL,
                    type='character',
                    help='A dictionary of arguments to pass on to the metric, such as the p value for Minkowski distance. If NULL then no arguments are passed on.')
,
    make_option(
                    c('--angular.rp.forest'),
                    action='store_true',
                    default=FALSE,
                    type='logical',
                    help='Whether to use an angular random projection forest to initialise the approximate nearest neighbor search. This can be faster, but is mostly on useful for metric that use an angular style distance such as cosine, correlation etc. In the case of those metrics angular forests will be chosen automatically.')
,
    make_option(
                    c('--do-not-verbose'),
                    action='store_false',
                    default=TRUE,
                    type='logical',
                    help='Controls verbosity')
,
    make_option(
                    c('--reduction.name'),
                    action='store',
                    default='umap',
                    type='character',
                    help='Name to store dimensional reduction under in the Seurat object')
,
    make_option(
                    c('--reduction.key'),
                    action='store',
                    default='UMAP',
                    type='character',
                    help='dimensional reduction key, specifies the string before the number for the dimension names. UMAP by default')
,
    make_option(
                    c('-o','--output-path'),
                    action='store',
                    type='character',
                    help='FILE IN')
,
    make_option(
                    c('--output-format'),
                    action='store',
                    default='seurat',
                    type='character',
                    help='FILE IN')
)

opt <- wsc_parse_args(option_list,
                      mandatory = c('input_path','output_path'))


if ( ! file.exists(opt$input_path) ) {
    stop((paste('File', opt$input_path, 'does not exist')))
}



if (! is.null(opt$dims) ) {
    opt$dims <- eval(parse(text=opt$dims))
}



features <- opt$features
if (! is.null(features) ) {
    features <- unlist(strsplit(opt$features, sep=','))
}



if (! is.null(opt$metric.kwds) ) {
    opt$metric.kwds <- eval(parse(text=opt$metric.kwds))
}


seurat_object <- read_seurat4_object(input_path = opt$input_path,
                    format = opt$input_format)

seurat_object_umap <- RunUMAP(object = seurat_object,
                    dims = opt$dims,
                    reduction = opt$reduction,
                    features = features,
                    graph = opt$graph,
                    assay = opt$assay,
                    nn.name = opt$nn.name,
                    slot = opt$slot,
                    umap.method = opt$umap.method,
                    reduction.model = opt$reduction.model,
                    return.model = opt$return.model,
                    n.neighbors = opt$n.neighbors,
                    n.components = opt$n.components,
                    metric = opt$metric,
                    n.epochs = opt$n.epochs,
                    learning.rate = opt$learning.rate,
                    min.dist = opt$min.dist,
                    spread = opt$spread,
                    set.op.mix.ratio = opt$set.op.mix.ratio,
                    local.connectivity = opt$local.connectivity,
                    repulsion.strength = opt$repulsion.strength,
                    negative.sample.rate = opt$negative.sample.rate,
                    a = opt$a,
                    b = opt$b,
                    uwot.sgd = opt$uwot.sgd,
                    seed.use = opt$seed.use,
                    metric.kwds = opt$metric.kwds,
                    angular.rp.forest = opt$angular.rp.forest,
                    verbose = !opt$do_not_verbose,
                    reduction.name = opt$reduction.name,
                    reduction.key = opt$reduction.key)

write_seurat4_object(seurat_object = seurat_object_umap,
                    output_path = opt$output_path,
                    format = opt$output_format)
