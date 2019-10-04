#!/usr/bin/env Rscript 

# This script currently calls an embedded function which is part of Seurat 3. As the package moves with Seurat 3, we will
# remove that function and rely on the native one.

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
    c("-o", "--output-directory"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'Seurat'.'"
  ),
  make_option(
    c("-n", "--study-name"),
    action= "store",
    default = "Seurat_study",
    type="character",
    help="Study name to be displayed in CellBrowser"
  ),
  make_option(
    c("-m", "--markers-file"),
    default = NULL,
    type="character",
    help="Path to markers file computed by Seurat"
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_directory'))
opt$study_name = gsub(" ","_", opt$study_name) # study names cannot have spaces.

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

suppressPackageStartupMessages(require(Seurat))
if(opt$input_format == "loom" ) {
  suppressPackageStartupMessages(require(loomR))
} else if(opt$input_format == "singlecellexperiment" ) {
  suppressPackageStartupMessages(require(scater))
}

# Input from serialized R object
seurat_object <- read_seurat3_object(input_path = opt$input_object_file, format = opt$input_format)

skip_markers = TRUE
if(!is.null(opt$markers_file)) {
  skip_markers = FALSE
}

ExportToCellbrowser(object = seurat_object, 
                    dir = opt$output_directory, 
                    markers.file = opt$markers_file, 
                    dataset.name = opt$study_name)
