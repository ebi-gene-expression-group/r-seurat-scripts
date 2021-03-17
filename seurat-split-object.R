#!/usr/bin/env Rscript

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))

# Load common functions

suppressPackageStartupMessages(require(workflowscriptscommon))

# parse options

option_list = list(
  make_option(
    c("-i", "--input-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = "RDS/Loom/SCE serialised object with content to split."
  ),
  make_option(
    c("-m", "--metadata-rds"),
    action = "store",
    default = NULL,
    type = 'character',
    help = 'Optional RDS with metadata.'
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
	      c("-o", "--output-path"),
	      action = "store",
	      default = "./",
	      type = "character",
	      help = "Output path, where split files will be left"
	      ),
  make_option(
    c("--split-by"),
    action = "store",
    type = "character",
    help = "Field to split data by. Required."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object', 'split_by'))


# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))
if(opt$input_format == "loom" | opt$output_format == "loom") {
  suppressPackageStartupMessages(require(loomR))
} else if(opt$input_format == "singlecellexperiment" | opt$output_format == "singlecellexperiment") {
  suppressPackageStartupMessages(require(scater))
}

seurat_object <- CreateSeuratObject(read_seurat3_object(input_path = opt$input_object, format = opt$input_format))
# Maybe check if the read object is not a Seurat object and call CreateSeuratObject if not.

if(!is.null(opt$metadata_rds)) {
  metadata<-readRDS(file=opt$metadata_rds)
  seurat_object <- AddMetaData(seurat_object, metadata=metadata)
}

output.list <- SplitObject(seurat_object, split.by = opt$split_by)

for(output in output.list) {
  write_seurat3_object(seurat_object = output,
                     output_path = paste0(opt$output_path, opt$split_by, output[opt$split_by][0], sep="_"),
                     format = opt$output_format)
}

