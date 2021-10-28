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
    c("-s", "--selection-method"),
    action = "store",
    default = "vst",
    type='character',
    help="How to choose top variable features. Choose one of: 'vst', 'mvp', disp."
  ),
  make_option(
    c("-m", "--mean-function"),
    action = "store",
    default = 'FastExpMean',
    type = 'character',
    help = "Function to compute x-axis value (average expression). Default is to take the mean of the detected (i.e. non-zero) values."
  ),
  make_option(
    c("-d", "--dispersion-function"),
    action = "store",
    default = 'FastLogVMR',
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
  ),
  make_option(
    c("--loess-span"),
    action = "store",
    default = 0.3,
    metavar = "Loess span parameter for VST",
    type = 'double',
    help = "(vst method) Loess span parameter used when fitting the variance-mean relationship. Default: 0.3"
  ),
  make_option(
    c("--clip-max"),
    action = "store",
    default = "auto",
    metavar = "Clip max for VST",
    type = 'character',
    help = "(vst method) After standardization values larger than clip.max will be set to clip.max; default is 'auto' which sets this value to the square root of the number of cells."
  ),
  make_option(
    c("--num-bin"),
    action = "store",
    default = 20,
    metavar = "Bins for scaled analysis",
    type = 'integer',
    help = "Total number of bins to use in the scaled analysis (default is 20)."
  ),
  make_option(
    c("--binning-method"),
    action = "store",
    default = "equal_width",
    metavar = "Binning method",
    type = 'character',
    help = "Specifies how the bins should be computed. Available methods are either equal_width: each bin is of equal width along the x-axis [default]; or equal_frequency: each bin contains an equal number of features (can increase statistical power to detect overdispersed features at high expression values, at the cost of reduced resolution along the x-axis)"
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file', 'output_text_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
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
# clean previous find variable genes
seurat_object@assays$RNA@meta.features<-data.frame(row.names = rownames(seurat_object@assays$RNA@meta.features))
variable_genes_seurat_object <- FindVariableFeatures(seurat_object, 
                                                  selection.method = opt$selection_method,
                                                  nfeatures = opt$nfeatures,
                                                  mean.function = opt$mean_function, 
                                                  dispersion.function = opt$dispersion_function, 
                                                  mean.cutoff = c(opt$x_low_cutoff, opt$x_high_cutoff),
                                                  dispersion.cutoff = c(opt$y_low_cutoff, opt$y_high_cutoff),
                                                  loess.span = opt$loess_span,
                                                  clip.max = opt$clip_max,
                                                  num.bin = opt$num_bin,
                                                  binning.method = opt$binning_method, 
                                                  verbose = FALSE)

# Output to a serialized R object
write_seurat4_object(seurat_object = variable_genes_seurat_object, 
                     output_path = opt$output_object_file, 
                     format = opt$output_format)

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
