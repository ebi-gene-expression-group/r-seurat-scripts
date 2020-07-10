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
    help = "Dimension for x-axis (default 1)"
  ),
  make_option(
    c("-b", "--dim-2"),
    action = "store",
    default = 2,
    type = 'integer',
    help = "Dimension for y-axis (default 2)"
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
    help = "Adjust point size for plotting"
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
    c("-m", "--do-bare"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "Do only minimal formatting (default : FALSE)"
  ),
  make_option(
    c("-u", "--cols-use"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Comma-separated list of colors, each color corresponds to an identity class. By default, ggplot assigns colors."
  ),
  make_option(
    c("-e", "--pt-shape"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "If NULL, all points are circles (default). You can specify any cell attribute (that can be pulled with FetchData) allowing for both different colors and different shapes on cells."
  ),
  make_option(
    c("-x", "--coord-fixed"),
    action = "store",
    default = FALSE,
    type = 'character',
    help = "Use a fixed scale coordinate system (for spatial coordinates). Default is FALSE"
  ),
  make_option(
    c("-n", "--no-axes"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "Setting to TRUE will remove the axes."
  ),              
  make_option(
    c("-k", "--dark-theme"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "Use a dark theme for the plot."
  ),              
  make_option(
    c("-q", "--plot-order"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Comma-separated string specifying the order of plotting for the idents (clusters). This can be useful for crowded plots if points of interest are being buried. Provide either a full list of valid idents or a subset to be plotted last (on top).."
  ),              
  make_option(
    c("-w", "--png-width"),
    action = "store",
    default = 1000,
    type = 'integer',
    help = "Width of png (px)."
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
if(opt$input_format == "loom" ) {
  suppressPackageStartupMessages(require(loomR))
} else if(opt$input_format == "singlecellexperiment" ) {
  suppressPackageStartupMessages(require(scater))
}

# Input from serialized R object

seurat_object <- read_seurat3_object(input_path = opt$input_object_file, format = opt$input_format)

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

# Parse color list (if present)

cols_use <- opt$cols_use
if (! is.null(cols_use)){
  cols_use <- wsc_split_string(cols_use)
}

# Parse plot order ident list (if present)

plot_order <- opt$plot_order
if (! is.null(plot_order)){
  plot_order <- wsc_split_string(plot_order)
}

# Open the image

png(filename = opt$output_image_file, width = opt$png_width, height = opt$png_height)
DimPlot(seurat_object, reduction.use = opt$reduction_use, dims = c(opt$dim_1, opt$dim_2), cells.use = cells_use, pt.size = opt$pt_size, label.size = opt$label_size, do.label = opt$do_label, group.by = opt$group_by, do.bare=opt$do_bare, cols.use=cols_use, pt.shape=opt$pt_shape, coord.fixed=opt$coord_fixed, no.axes=opt$no_axes, dark.theme=opt$dark_theme, plot.order=plot_order, plot.title = opt$plot_title)
dev.off()
