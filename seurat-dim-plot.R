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
    c("-r", "--reduction-use"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Which dimensionality reduction to use. Default is "pca", can also be "tsne", or "ica", assuming these are precomputed.'
  ),
  make_option(
    c("-a", "--dim-1"),
    action = "store",
    default = 1,
    type = 'integer',
    help = "File name in which a serialized R matrix object may be found."
  ),
  make_option(
    c("-b", "--dim-2"),
    action = "store",
    default = 2,
    type = 'integer',
    help = "File name in which a serialized R matrix object may be found."
  ),
  make_option(
    c("-c", "--cells-use"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "File to be used to derive a vector of cells to plot (default is all cells)."
  ),
  make_option(
    c("-p", "--pt-size"),
    action = "store",
    default = 1,
    type = 'integer',
    help = "File to be used to derive a vector of cells to plot (default is all cells)."
  ),
  make_option(
    c("-l", "--label-size"),
    action = "store",
    default = 4,
    type = 'integer',
    help = "Sets size of labels."
  ),
  make_option(
    c("-d", "--do-label"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "Whether to label the clusters."
  ),
  make_option(
    c("-f", "--group-by"),
    action = "store",
    default = 'ident',
    type = 'character',
    help = "Group (color) cells in different ways (for example, orig.ident)."
  ),
  make_option(
    c("-t", "--plot-title"),
    action = "store",
    default = 1,
    type = 'character',
    help = "Title for plot."
  ),
  make_option(
    c("-w", "--png-width"),
    action = "store",
    default = 1000,
    type = 'integer',
    help = "Widthof png (px)."
  ),
  make_option(
    c("-j", "--png-height"),
    action = "store",
    default = 1000,
    type = 'integer',
    help = "Height of png file (px)."
  ),
  make_option(
    c("-o", "--output-image-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to save the PCA image"
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_image_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))

# Input from serialized R object

seurat_object <- readRDS(opt$input_object_file)

# Read cells file (if present)

if (! is.null(opt$cells_use)){
  if (! file.exists(opt$cells_use)){
    stop((paste('Supplied genes file', opt$cells_use, 'does not exist')))
  }else{
    cells_use <- readLines(opt$cells_use)
  }
}else{
  cells_use <- NULL
}

# Open the image

png(filename = opt$output_image_file, width = opt$png_width, height = opt$png_height)
DimPlot(seurat_object, reduction.use = opt$reduction_use, dim.1 = opt$dim_1, dim.2 = opt$dim_2, cells.use = cells_use, pt.size = opt$pt_size, label.size = opt$label_size, do.label = opt$do_label, group.by = opt$group_by, plot.title = opt$plot_title)
dev.off()