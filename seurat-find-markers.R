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
    help = "File to be used to derive a vector of genes to test. Default is to use all genes."
  ),
  make_option(
    c("--input-format"),
    action = "store",
    default = "seurat",
    type = 'character',
    help = "Either loom, seurat, anndata or singlecellexperiment for the input format to read."
  ),
  make_option(
    c("-l", "--logfc-threshold"),
    action = "store",
    default = 0.25,
    type = 'double',
    help = "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing logfc.threshold speeds up the function, but can miss weaker signals."
  ),
  make_option(
    c("-m", "--min-pct"),
    action = "store",
    default = 0.1,
    type = 'double',
    help = "Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1."
  ),
  make_option(
    c("-n", "--min-diff-pct"),
    action = "store",
    default = -Inf,
    type = 'double',
    help = "Only test genes that show a minimum difference in the fraction of detection between the two groups. Set to -Inf by default."
  ),
  make_option(
    c("-p", "--only-pos"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "Only return positive markers (FALSE by default)."
  ),
  make_option(
    c("-t", "--test-use"),
    action = "store",
    default = 'wilcox',
    type = 'character',
    help = "Denotes which test to use. Available options are 'wilcox', 'bimod', 'roc', 't', 'tobit', 'poisson', 'negbinom', 'MAST', 'DESeq2'. See ?FindMarkers() for more info."
  ),
  make_option(
    c("-x", "--max-cells-per-ident"),
    action = "store",
    default = Inf,
    type = 'double',
    help = "Down sample each identity class to a max number. Default is no downsampling. Not activated by default (set to Inf)."
  ),
  make_option(
    c("-c", "--min-cells-gene"),
    action = "store",
    default = 3,
    type = 'integer',
    help = "Minimum number of cells expressing the gene in at least one of the two groups, currently only used for poisson and negative binomial tests."
  ),
  make_option(
    c("-d", "--min-cells-group"),
    action = "store",
    default = 3,
    type = 'integer',
    help = "Minimum number of cells in one of the groups."
  ),
  make_option(
    c("-o", "--output-text-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store text format matrix containing a ranked list of putative markers, and associated statistics (p-values, ROC score, etc.)."
  )
  # , # The defaults for the above variables lead to zero marker genes found:
  # make_option(
  #   c("--slot"),
  #   action = "store",
  #   default = NULL,
  #   metavar = "Data slot",
  #   type = 'character',
  #   help = "Slot to pull data from; note that if test.use is 'negbinom', 'poisson', or 'DESeq2', slot will be set to 'counts'"
  # ),
  # make_option(
  #   c("--latent-vars"),
  #   action = "store",
  #   default = NULL,
  #   metavar = "Data slot",
  #   type = 'character',
  #   help = "Comma separated variables to test, used only when test.use is one of 'LR', 'negbinom', 'poisson', or 'MAST'"
  # ),
  # make_option(
  #   c("--pseudocount-use"),
  #   action = "store",
  #   default = NULL,
  #   metavar = "Use pseudocount",
  #   type = 'integer',
  #   help = "Pseudocount to add to averaged expression values when calculating logFC. 1 by default."
  # )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_text_file'))

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

# Latent vars
#latent_vars<-NULL
#if (! is.null(opt$latent_vars)) {
#  latent_vars<-strsplit(opt$latent_vars,split = ",")
#}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))
if(opt$input_format == "loom" ) {
  suppressPackageStartupMessages(require(SeuratDisk))
} else if(opt$input_format == "singlecellexperiment" ) {
  suppressPackageStartupMessages(require(scater))
}

# Input from serialized R object

seurat_object <- read_seurat4_object(input_path = opt$input_object_file, format = opt$input_format)

# Get results matrix
results_matrix<-FindAllMarkers(
  seurat_object, 
  # slot = opt$slot,
  features = genes_use, 
  # latent.vars = latent_vars, 
  # pseudocount.use = opt$pseudocount_use,
  logfc.threshold = opt$logfc_threshold,
  test.use = opt$test_use,
  min.pct = opt$min_pct,
  min.diff.pct = opt$min_diff_pct,
  verbose = FALSE,
  min.cells.feature = opt$min_cells_gene,
  min.cells.group = opt$min_cells_group,
  only.pos = opt$only_pos,
  max.cells.per.ident = opt$max_cells_per_ident
)

# Summarise output

# Some parameters aren't interesting for reporting purposes (e.g. file
# locations), so hide from the summary

nonreport_params <- c('input_object_file', 'output_object_file', 'help', 'output_text_file')
opt_table <- data.frame(value=unlist(opt), stringsAsFactors = FALSE)
opt_table <- opt_table[! rownames(opt_table) %in% nonreport_params, , drop = FALSE]

markers_by_cluster <- merge(
  data.frame(as.matrix(table(Idents(seurat_object)))),
  data.frame(as.matrix(table(results_matrix$cluster))),
  by = 'row.names'
)
rownames(markers_by_cluster) <- markers_by_cluster[,1]
colnames(markers_by_cluster) <- c('Cluster', 'No. cells', 'No. markers')

cat(c(
  "Summary of markers found:\n",
  capture.output(markers_by_cluster[,-1]),
  '\nParameter values:', 
  capture.output(print(opt_table))
  ), sep = '\n')

# Output variable genes to a simple text file

write.csv(results_matrix, file = opt$output_text_file, row.names = TRUE)
