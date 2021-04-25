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
    help = "Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 Leiden)."
  ),
  make_option(
    c("-m", "--tmp-file-location"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Directory where intermediate files will be written. Specify the ABSOLUTE path."
  ),
  make_option(
    c("--modularity-fxn"),
    action = "store",
    default = 1,
    type = 'integer',
    help = "Modularity function: 1 standard, 2 alternative."
  ),
  make_option(
    c("--method"),
    action = "store",
    default = 'matrix',
    type = 'character',
    help = "Method for leiden  (defaults to matrix which is fast for small datasets). Enable method = \"igraph\" to avoid casting large data to a dense matrix."
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
  ),
  make_option(
    c("--graph-name"),
    action = "store",
    default = "RNA_nn",
    type = 'character',
    help = "Name of graph to use for the clustering algorithm."
  ),
  make_option(
    c("-s", "--nrandom-starts"),
    action = "store",
    default = 10,
    type = 'integer',
    help = "Number of random starts"
  ),
  make_option(
    c("--n-iterations"),
    action = "store",
    default = 10,
    type = 'integer',
    help = "Maximal number of iterations per random start"
  ),
  make_option(
    c("--group-singletons"),
    action = "store_true",
    default = FALSE,
    help = "Group singletons into nearest cluster. If FALSE, assign all singletons to a \"singleton\" group"
  ),
  make_option(
    c("--random-seed"),
    action = "store",
    default = 0,
    type = 'integer',
    help = "Seed of the random number generator"
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file', 'output_text_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
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
  suppressPackageStartupMessages(require(SeuratDisk))
} else if(opt$input_format == "singlecellexperiment" | opt$output_format == "singlecellexperiment") {
  suppressPackageStartupMessages(require(scater))
}

# Input from serialized R object

seurat_object <- read_seurat4_object(input_path = opt$input_object_file, format = opt$input_format)

clustered_object <- FindClusters(seurat_object, 
                                 algorithm = opt$algorithm,
                                 modularity.fxn = opt$modularity_fxn,
                                 method = opt$method,
                                 n.start = opt$nrandom_starts,
                                 n.iter = opt$n_iterations,
                                 random.seed = opt$random_seed,
                                 group.singletons = opt$group_singletons,
                                 verbose = TRUE, 
                                 resolution = opt$resolution, 
                                 graph.name = opt$graph_name,
                                 temp.file.location = opt$temp_file_location)

# Summarise the clustering

# Some parameters aren't interesting for reporting purposes (e.g. file
# locations), so hide from the summary

nonreport_params <- c('input_object_file', 'output_object_file', 'help', 'output_text_file', 'tmp_file_location')
opt_table <- data.frame(value=unlist(opt), stringsAsFactors = FALSE)
opt_table <- opt_table[! rownames(opt_table) %in% nonreport_params, , drop = FALSE]

cluster_table <- as.data.frame(table(Idents(clustered_object)))
colnames(cluster_table) <- c('Cluster', 'No. cells')
rownames(cluster_table) <- cluster_table$Cluster

cat(paste(ncol(GetAssayData(clustered_object)), 'cells fall into ', length(levels(Idents(clustered_object))), 'final clusters. Membership numbers:\n'), capture.output(cluster_table[,2, drop = FALSE]), '\nParameter values:\n', capture.output(print(opt_table)), sep = '\n')

# Output to a serialized R object
write_seurat4_object(seurat_object = clustered_object, 
                     output_path = opt$output_object_file, 
                     format = opt$output_format)

# Output variable genes to a simple text file

write.csv(data.frame(cell=names(Idents(clustered_object)), cluster=Idents(clustered_object)), file = opt$output_text_file, row.names = FALSE)
