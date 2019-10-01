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
    c("-s", "--selection-method"),
    action = "store",
    default = "vst",
    type='character',
    help="How to choose top variable features. Choose one of: 'vst', 'mvp', disp."
  ),
  make_option(
    c("-m", "--mean-function"),
    action = "store",
    default = 'ExpMean',
    type = 'character',
    help = "Function to compute x-axis value (average expression). Default is to take the mean of the detected (i.e. non-zero) values."
  ),
  make_option(
    c("-d", "--dispersion-function"),
    action = "store",
    default = 'LogVMR',
    type = 'character',
    help = "Function to compute y-axis value (dispersion). Default is to take the standard deviation of all values."
  ),
  make_option(
    c("-l", "--x-low-cutoff"),
    action = "store",
    default = 0.1,
    type = 'double',
    help = "Bottom cutoff on x-axis (mean) for identifying variable genes."
  ),
  make_option(
    c("-j", "--x-high-cutoff"),
    action = "store",
    default = 8,
    type = 'double',
    help = "Top cutoff on x-axis (mean) for identifying variable genes."
  ),
  make_option(
    c("-n", "--nfeatures"),
    action = "store",
    default = 2000,
    type = 'integer',
    help = "Number of features to return."
  ),
  make_option(
    c("-y", "--y-low-cutoff"),
    action = "store",
    default = 1,
    type = 'double',
    help = "Bottom cutoff on y-axis (dispersion) for identifying variable genes."
  ),
  make_option(
    c("-z", "--y-high-cutoff"),
    action = "store",
    default = Inf,
    type = 'double',
    help = "Top cutoff on y-axis (dispersion) for identifying variable genes."
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
    help = "File name in which to store variable genes in plain text."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file', 'output_text_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))

# Input from serialized R object

seurat_object <- readRDS(opt$input_object_file)
variable_genes_seurat_object <- FindVariableFeatures(seurat_object, 
                                                  selection.method = opt$selection_method,
                                                  nfeatures = opt$nfeatures,
                                                  mean.function = get(opt$mean_function), 
                                                  dispersion.function = get(opt$dispersion_function), 
                                                  mean.cutoff = c(opt$x_low_cutoff, opt$x_high_cutoff),
                                                  dispersion.cutoff = c(opt$y_low_cutoff, opt$y_high_cutoff),
                                                  verbose = FALSE)

# Output to a serialized R object

saveRDS(variable_genes_seurat_object, file = opt$output_object_file)

# Output variable gene numbers

# Some parameters aren't interesting for reporting purposes (e.g. file
# locations), so hide from the summary

nonreport_params <- c('input_object_file', 'output_object_file', 'help', 'output_text_file')
opt_table <- data.frame(value=unlist(opt), stringsAsFactors = FALSE)
opt_table <- opt_table[! rownames(opt_table) %in% nonreport_params, , drop = FALSE]

cat(c(
    paste(length(VariableFeatures(variable_genes_seurat_object)), 'variable genes detected out of total', nrow(GetAssayData(seurat_object)))),
    '\nParameter values:', 
    capture.output(print(opt_table)
), sep = '\n')          

# Output variable genes to a simple text file

writeLines(con=opt$output_text_file, VariableFeatures(variable_genes_seurat_object))
