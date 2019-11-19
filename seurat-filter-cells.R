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
    c("-s", "--subset-names"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Parameters to subset on. Eg, the name of a gene, PC1, a column name in object@meta.data, etc. Any argument that can be retreived using FetchData."
  ),
  make_option(
    c("-l", "--low-thresholds"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Low cutoffs for the parameters (default is -Inf)."
  ),
  make_option(
    c("-j", "--high-thresholds"),
    action = "store",
    default = NA,
    type = 'character',
    help = "High cutoffs for the parameters (default is Inf)."
  ),
  make_option(
    c("-c", "--cells-use"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Comma-separated list of cell names to use as a subset. Alternatively, text file with one cell per line."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'Seurat'.'"
  ),
  make_option(
    c("--idents"),
    action = "store",
    default = NULL,
    type = 'character',
    metavar = 'Ident classes to keep',
    help = "Comma-separated list of identity classes to keep"
  ),
  make_option(
    c("--features"),
    action = "store",
    default = NULL,
    type = 'character',
    metavar = 'Features to keep',
    help = "Comma-separated list or file path with features (normally genes) to keep"
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'subset_names', 'output_object_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

# Parse numeric fields
subset_names <- wsc_split_string(opt$subset_names)
lt <- wsc_parse_numeric(opt, 'low_thresholds', -Inf, length(subset_names))
ht <- wsc_parse_numeric(opt, 'high_thresholds', Inf, length(subset_names))

# Check the cells_use

cells_use <- opt$cells_use
if (! is.null(cells_use)){
  if (file.exists(cells_use)){
    cells_use <- read.table(file = cells_use, stringsAsFactors = FALSE)$V1
  }else{
    cells_use <- wsc_split_string(cells_use)
  }
  
  check_cells(cells_use, seurat_object)
  print(paste('Filtering to', length(cells_use), 'specified cells'))
}

# Check features to use
features_use<-NULL
if (! is.null(opt$features) ) {
  if (file.exists(opt$features)){
    features_use <- read.table(file = opt$features, stringsAsFactors = FALSE)$V1
  }else{
    features_use <- wsc_split_string(opt$features)
  }
}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))
if(opt$input_format == "loom" | opt$output_format == "loom") {
  suppressPackageStartupMessages(require(loomR))
} else if(opt$input_format == "singlecellexperiment" | opt$output_format == "singlecellexperiment") {
  suppressPackageStartupMessages(require(scater))
}

# Input from serialized R object

seurat_object <- read_seurat3_object(input_path = opt$input_object_file, format = opt$input_format)
# Are the metadata variables valid for this object?
wsc_check_metadata(seurat_object, subset_names)

# Given the new setup, now we need to iterate over all the elements provided for filtering
# and come up with an intersection on all of them
cells_bool<-rep(TRUE, length(subset_names))
for(i in seq_along(subset_names)) {
  cells_bool<-cells_bool & seurat_object[[]][subset_names[[i]]] > lt[i] & seurat_object[[]][subset_names[[i]]] < ht[i]
  cat(c(length(which(cells_bool)),"cells remaining after applying",subset_names[i],"thresholds",lt[i],"< x <",ht[i],"\n"))
}
# Intersect that with any explicit cells that the user has asked for, if any
cells<-colnames(seurat_object)[cells_bool]
if (! is.null(cells_use)){
  cells<-intersect(cells_use, cells)
  cat(c(lenght(cells)," cells remaining after intersecting threshold based selection with provided list of cells.\n"))
}

filtered_seurat_object <- subset(seurat_object, cells = cells)

if (! is.null(features_use) ) {
  filtered_seurat_object <- subset(filtered_seurat_object, features=features_use)
}
if (! is.null(opt$idents) ) {
  idents_vector=unlist(strsplit(opt$idents, split=","))
  filtered_seurat_object <- subset(filtered_seurat_object, idents = idents_vector)
}

# Print a summary of the affects of filtering

# Some parameters aren't interesting for reporting purposes (e.g. file
# locations), so hide from the summary

nonreport_params <- c('input_object_file', 'output_object_file', 'help')
opt_table <- data.frame(value=unlist(opt), stringsAsFactors = FALSE)
opt_table <- opt_table[! rownames(opt_table) %in% nonreport_params, , drop = FALSE]

cat(c(
  '# Before filtering:', 
  capture.output(print(seurat_object)), 
  '\n# After filtering:', 
  capture.output(print(filtered_seurat_object)),
  '\n# Parameter values:\n',
  capture.output(print(opt_table))
), 
sep = '\n')

# Output to a serialized R object
write_seurat3_object(seurat_object = filtered_seurat_object, 
                     output_path = opt$output_object_file, 
                     format = opt$output_format)
