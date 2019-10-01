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
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'subset_names', 'output_object_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

# Input from serialized R object

seurat_object <- readRDS(opt$input_object_file)

# Are the metadata variables valid for this object?

subset_names <- wsc_split_string(opt$subset_names)
wsc_check_metadata(seurat_object, subset_names)

# Parse numeric fields

lt <- wsc_parse_numeric(opt, 'low_thresholds', -Inf, length(subset_names))
ht <- wsc_parse_numeric(opt, 'high_thresholds', Inf, length(subset_names))

# Check the cells_use

cells_use <- opt$cells_use
if (! is.null(cells_use)){
  if (file.exists(cells_use)){
    cells_use <- read.table("test_cells.txt", stringsAsFactors = FALSE)$V1
  }else{
    cells_use <- wsc_split_string(cells_use)
  }
  
  check_cells(cells_use, seurat_object)
  print(paste('Filtering to', length(cells_use), 'specified cells'))
}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))

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

saveRDS(filtered_seurat_object, file = opt$output_object_file)
