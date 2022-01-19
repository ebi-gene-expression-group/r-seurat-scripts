#!/usr/bin/env Rscript
# This script has been automatically generated through
#
# YAML2RScript.py -i ../r-seurat-scripts/cli-generated/manually_crafted_YAML/seurat-plot.yaml -o ../r-seurat-scripts/seurat-plot.R
#
# to change this file edit the input YAML and re-run the above command

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(SeuratDisk))
suppressPackageStartupMessages(require(patchwork))

option_list <- list(
    make_option(
        c("-i", "--input-object-file"),
        action = "store",
        metavar = "Input file",
        type = "character",
        help = "Query file with Seurat object in either RDS-Seurat, Loom or SCE"
    ),
    make_option(
        c("--input-format"),
        action = "store",
        default = "seurat",
        metavar = "Input format",
        type = "character",
        help = "Either loom, seurat, anndata or singlecellexperiment for the input format to read."
    ),
    make_option(
        c("--plot-type"),
        action = "store",
        default = "FeaturePlot",
        type = "character",
        help = "Either FeaturePlot, RidgePlot, DimPlot, VlnPlot or DotPlot."
    ),
    make_option(
        c("--dims"),
        action = "store",
        default = "1,2",
        type = "character",
        help = "Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions"
    ),
    make_option(
        c("--cells"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Vector of cells to plot (default is all cells)"
    ),
    make_option(
        c("--cols"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Vector of colors, each color corresponds to an identity class. This may also be a single character or numeric value corresponding to a palette as specified by brewer.pal.info. By default, ggplot2 assigns colors. We also include a number of palettes from the pals package. See 'DiscretePalette' for details."
    ),
    make_option(
        c("--pt-size"),
        action = "store",
        default = NULL,
        type = "integer",
        help = "Adjust point size for plotting"
    ),
    make_option(
        c("--reduction"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca"
    ),
    make_option(
        c("--group-by"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Name of one or more metadata columns to group (color) cells by (for example, orig.ident); pass 'ident' to group by identity class"
    ),
    make_option(
        c("--split-by"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Name of a metadata column to split plot by; see 'FetchData' for more details"
    ),
    make_option(
        c("--shape-by"),
        action = "store",
        default = NULL,
        type = "character",
        help = "If NULL, all points are circles (default). You can specify any cell attribute (that can be pulled with FetchData) allowing for both different colors and different shapes on cells.  Only applicable if raster is FALSE."
    ),
    make_option(
        c("--order"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Specify the order of plotting for the idents. This can be useful for crowded plots if points of interest are being buried. Provide either a full list of valid idents or a subset to be plotted last (on top)"
    ),
    make_option(
        c("--shuffle"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Whether to randomly shuffle the order of points. This can be useful for crowded plots if points of interest are being buried. (default is FALSE)"
    ),
    make_option(
        c("--seed"),
        action = "store",
        default = 1,
        type = "integer",
        help = "Sets the seed if randomly shuffling the order of points."
    ),
    make_option(
        c("--label"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Whether to label the clusters"
    ),
    make_option(
        c("--label-size"),
        action = "store",
        default = 4,
        type = "integer",
        help = "Sets size of labels"
    ),
    make_option(
        c("--label-color"),
        action = "store",
        default = "black",
        type = "character",
        help = "Sets the color of the label text"
    ),
    make_option(
        c("--label-box"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Whether to put a box around the label text (geom_text vs geom_label)"
    ),
    make_option(
        c("--repel"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Repel labels"
    ),
    make_option(
        c("--cells-highlight"),
        action = "store",
        default = NULL,
        type = "character",
        help = "A list of character or numeric vectors of cells to highlight. If only one group of cells desired, can simply pass a vector instead of a list. If set, colors selected cells to the color(s) in  'cols.highlight'  and other cells black (white if dark.theme = TRUE); will also resize to the size(s) passed to 'sizes.highlight'."
    ),
    make_option(
        c("--cols-highlight"),
        action = "store",
        default = "#DE2D26",
        type = "character",
        help = "A vector of colors to highlight the cells as; will repeat to the length groups in cells.highlight. Comma separated."
    ),
    make_option(
        c("--sizes-highlight"),
        action = "store",
        default = 1,
        type = "integer",
        help = "Size of highlighted cells; will repeat to the length groups in cells.highlight"
    ),
    make_option(
        c("--na-value"),
        action = "store",
        default = "grey50",
        type = "character",
        help = "Color value for NA points when using custom scale"
    ),
    make_option(
        c("--ncol"),
        action = "store",
        default = NULL,
        type = "integer",
        help = "Number of columns for display when combining plots"
    ),
    make_option(
        c("--features"),
        action = "store",
        default = "",
        type = "character",
        help = "Features to plot (gene expression, metrics, PC scores, anything that can be retreived by FetchData)"
    ),
    make_option(
        c("--cols-ridgplot"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Colors to use for plotting, comma separated"
    ),
    make_option(
        c("--idents"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Which classes to include in the plot (default is all)"
    ),
    make_option(
        c("--sort"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Sort identity classes (on the x-axis) by the average expression of the attribute being potted, can also pass 'increasing' or 'decreasing' to change sort direction"
    ),
    make_option(
        c("--assay"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Name of assay to use, defaults to the active assay"
    ),
    make_option(
        c("--y-max"),
        action = "store",
        default = NULL,
        type = "double",
        help = "Maximum y axis value"
    ),
    make_option(
        c("--same-y-lims"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Set all the y-axis limits to the same values"
    ),
    make_option(
        c("--log"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "plot the feature axis on log scale"
    ),
    make_option(
        c("--slot"),
        action = "store",
        default = "data",
        type = "character",
        help = "Use non-normalized counts data for plotting"
    ),
    make_option(
        c("--stack"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Horizontally stack plots for each feature"
    ),
    make_option(
        c("--fill-by"),
        action = "store",
        default = "feature",
        type = "character",
        help = "Color violins/ridges based on either 'feature' or 'ident'"
    ),
    make_option(
        c("--cols-feature-plot"),
        action = "store",
        default = "lightgrey,blue",
        type = "character",
        help = "The two colors to form the gradient over. Provide as string vector with the first color corresponding to low values, the second to high. Also accepts a Brewer color scale or vector of colors."
    ),
    make_option(
        c("--min-cutoff"),
        action = "store",
        default = "NA",
        type = "character",
        help = "Vector of minimum cutoff values for each feature, may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')"
    ),
    make_option(
        c("--max-cutoff"),
        action = "store",
        default = "NA",
        type = "character",
        help = "Vector of maximum cutoff values for each feature, may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')"
    ),
    make_option(
        c("--keep-scale"),
        action = "store",
        default = "feature",
        type = "character",
        help = "How to handle the color scale across multiple plots."
    ),
    make_option(
        c("--blend"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Scale and blend expression values to visualize coexpression of two features"
    ),
    make_option(
        c("--blend-threshold"),
        action = "store",
        default = 0,
        type = "integer",
        help = "The color cutoff from weak signal to strong signal; ranges from 0 to 1."
    ),
    make_option(
        c("--coord-fixed"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Plot cartesian coordinates with fixed aspect ratio"
    ),
    make_option(
        c("--do-not-by-col"),
        action = "store_false",
        default = TRUE,
        type = "logical",
        help = "If splitting by a factor, plot the splits per column with the features as rows; ignored if blend = TRUE."
    ),
    make_option(
        c("--adjust"),
        action = "store",
        default = 1,
        type = "integer",
        help = "Adjust parameter for geom_violin"
    ),
    make_option(
        c("--split-plot"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "plot each group of the split violin plots by multiple or single violin shapes."
    ),
    make_option(
        c("--flip"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "flip plot orientation (identities on x-axis)"
    ),
    make_option(
        c("--cols-dot-plot"),
        action = "store",
        default = "yellow,lightgrey,blue",
        type = "character",
        help = "Colors to plot: the name of a palette from RColorBrewer::brewer.pal.info , a pair of colors defining a gradient, or 3+ colors defining multiple gradients (if split.by is set)"
    ),
    make_option(
        c("--col-min"),
        action = "store",
        default = -2,
        type = "integer",
        help = "Minimum scaled average expression threshold (everything smaller will be set to this)"
    ),
    make_option(
        c("--col-max"),
        action = "store",
        default = 2,
        type = "integer",
        help = "Maximum scaled average expression threshold (everything larger will be set to this)"
    ),
    make_option(
        c("--dot-min"),
        action = "store",
        default = 0,
        type = "integer",
        help = "The fraction of cells at which to draw the smallest dot (default is 0). All cell groups with less than this expressing the given gene will have no dot drawn."
    ),
    make_option(
        c("--dot-scale"),
        action = "store",
        default = 6,
        type = "integer",
        help = "Scale the size of the points, similar to cex"
    ),
    make_option(
        c("--cluster-idents"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Whether to order identities by hierarchical clusters based on given features, default is FALSE"
    ),
    make_option(
        c("--do-not-scale"),
        action = "store_false",
        default = TRUE,
        type = "logical",
        help = "Determine whether the data is scaled, will be done by default"
    ),
    make_option(
        c("--scale-by"),
        action = "store",
        default = "radius",
        type = "character",
        help = "Scale the size of the points by 'size' or by 'radius'"
    ),
    make_option(
        c("--scale-min"),
        action = "store",
        default = "NA",
        type = "double",
        help = "Set lower limit for scaling, use NA for default"
    ),
    make_option(
        c("--scale-max"),
        action = "store",
        default = "NA",
        type = "double",
        help = "Set upper limit for scaling, use NA for default"
    ),
    make_option(
        c("--output-rds-file"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Path to RDS to keep the plot object"
    ),
    make_option(
        c("--plot-out"),
        action = "store",
        type = "character",
        help = "File name to create on disk."
    ),
    make_option(
        c("--plot-format"),
        action = "store",
        default = "png",
        type = "character",
        help = "Format to use, either PNG, EPS, PostScript, TeX, PDF, JPEG, TIFF or SVG"
    ),
    make_option(
        c("--scale-factor"),
        action = "store",
        default = 1,
        type = "double",
        help = "Multiplicative scaling factor."
    ),
    make_option(
        c("--width"),
        action = "store",
        default = 20,
        type = "double",
        help = "Width of the figure, in the selected units."
    ),
    make_option(
        c("--height"),
        action = "store",
        default = 20,
        type = "double",
        help = "Height of the figure, in the selected units."
    ),
    make_option(
        c("--units"),
        action = "store",
        default = "cm",
        type = "character",
        help = "Units for the plot dimensions."
    ),
    make_option(
        c("--dpi"),
        action = "store",
        default = 300,
        type = "integer",
        help = "Plot resolution. Also accepts a string input: retina (320), print (300), or screen (72). Applies only to raster output types."
    ),
    make_option(
        c("--do-not-limitsize"),
        action = "store_false",
        default = TRUE,
        type = "logical",
        help = "When TRUE (the default) ggsave() will not save images larger than 50x50 inches, to prevent the common error of specifying dimensions in pixels."
    ),
    make_option(
        c("--bg"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Background colour. If NULL, uses the plot.background fill value from the plot theme."
    )
)

opt <- wsc_parse_args(option_list, 
                      mandatory = c("input_object_file", "plot_out"))
                

if (!file.exists(opt$input_object_file)) {
    stop((paste("File", opt$input_object_file, "does not exist")))
}



dims <- opt$dims
if (!is.null(dims)) {
    dims <- as.numeric(unlist(strsplit(opt$dims, split = ",")))
}



cells <- opt$cells
if (!is.null(cells)) {
    cells <- unlist(strsplit(opt$cells, split = ","))
}



cols <- opt$cols
if (!is.null(cols)) {
    cols <- unlist(strsplit(opt$cols, split = ","))
}



group_by <- opt$group_by
if (!is.null(group_by)) {
    group_by <- unlist(strsplit(opt$group_by, split = ","))
}



order <- opt$order
if (!is.null(order)) {
    order <- unlist(strsplit(opt$order, split = ","))
}



cells_highlight <- opt$cells_highlight
if (!is.null(cells_highlight)) {
    cells_highlight <- unlist(strsplit(opt$cells_highlight, split = ","))
}



cols_highlight <- opt$cols_highlight
if (!is.null(cols_highlight)) {
    cols_highlight <- unlist(strsplit(opt$cols_highlight, split = ","))
}



features <- opt$features
if (!is.null(features)) {
    features <- unlist(strsplit(opt$features, split = ","))
}



cols_ridgplot <- opt$cols_ridgplot
if (!is.null(cols_ridgplot)) {
    cols_ridgplot <- unlist(strsplit(opt$cols_ridgplot, split = ","))
}



idents <- opt$idents
if (!is.null(idents)) {
    idents <- unlist(strsplit(opt$idents, split = ","))
}



group_by <- opt$group_by
if (!is.null(group_by)) {
    group_by <- unlist(strsplit(opt$group_by, split = ","))
}



features <- opt$features
if (!is.null(features)) {
    features <- unlist(strsplit(opt$features, split = ","))
}



dims <- opt$dims
if (!is.null(dims)) {
    dims <- as.numeric(unlist(strsplit(opt$dims, split = ",")))
}



cells <- opt$cells
if (!is.null(cells)) {
    cells <- unlist(strsplit(opt$cells, split = ","))
}



cols_feature_plot <- opt$cols_feature_plot
if (!is.null(cols_feature_plot)) {
    cols_feature_plot <- unlist(strsplit(opt$cols_feature_plot, split = ","))
}



min_cutoff <- opt$min_cutoff
if (!is.null(min_cutoff)) {
    min_cutoff <- unlist(strsplit(opt$min_cutoff, split = ","))
}



max_cutoff <- opt$max_cutoff
if (!is.null(max_cutoff)) {
    max_cutoff <- unlist(strsplit(opt$max_cutoff, split = ","))
}



features <- opt$features
if (!is.null(features)) {
    features <- unlist(strsplit(opt$features, split = ","))
}



cols <- opt$cols
if (!is.null(cols)) {
    cols <- unlist(strsplit(opt$cols, split = ","))
}



group_by <- opt$group_by
if (!is.null(group_by)) {
    group_by <- unlist(strsplit(opt$group_by, split = ","))
}



features <- opt$features
if (!is.null(features)) {
    features <- unlist(strsplit(opt$features, split = ","))
}



cols_dot_plot <- opt$cols_dot_plot
if (!is.null(cols_dot_plot)) {
    cols_dot_plot <- unlist(strsplit(opt$cols_dot_plot, split = ","))
}



idents <- opt$idents
if (!is.null(idents)) {
    idents <- unlist(strsplit(opt$idents, split = ","))
}



group_by <- opt$group_by
if (!is.null(group_by)) {
    group_by <- unlist(strsplit(opt$group_by, split = ","))
}


load_seurat4_packages_for_format(formats = c(opt$query_format, opt$anchors_format, opt$reference_format))

seurat_object <- read_seurat4_object(input_path = opt$input_object_file,
                    format = opt$input_format)
if ( opt$plot_type == 'DimPlot' ) { 
plot_object <- DimPlot(object = seurat_object,
                    dims = dims,
                    cells = cells,
                    cols = cols,
                    pt.size = opt$pt_size,
                    reduction = opt$reduction,
                    group.by = group_by,
                    split.by = opt$split_by,
                    shape.by = opt$shape_by,
                    order = order,
                    shuffle = opt$shuffle,
                    seed = opt$seed,
                    label = opt$label,
                    label.size = opt$label_size,
                    label.color = opt$label_color,
                    label.box = opt$label_box,
                    repel = opt$repel,
                    cells.highlight = cells_highlight,
                    cols.highlight = cols_highlight,
                    sizes.highlight = opt$sizes_highlight,
                    na.value = opt$na_value,
                    ncol = opt$ncol)
} else if ( opt$plot_type == 'RidgePlot' ) { 
plot_object <- RidgePlot(object = seurat_object,
                    features = features,
                    cols = cols_ridgplot,
                    idents = idents,
                    sort = opt$sort,
                    assay = opt$assay,
                    group.by = group_by,
                    y.max = opt$y_max,
                    same.y.lims = opt$same_y_lims,
                    log = opt$log,
                    ncol = opt$ncol,
                    slot = opt$slot,
                    stack = opt$stack,
                    fill.by = opt$fill_by)
} else if ( opt$plot_type == 'FeaturePlot' ) { 
plot_object <- FeaturePlot(object = seurat_object,
                    features = features,
                    dims = dims,
                    cells = cells,
                    cols = cols_feature_plot,
                    pt.size = opt$pt_size,
                    order = opt$order,
                    min.cutoff = min_cutoff,
                    max.cutoff = max_cutoff,
                    reduction = opt$reduction,
                    split.by = opt$split_by,
                    keep.scale = opt$keep_scale,
                    shape.by = opt$shape_by,
                    slot = opt$slot,
                    blend = opt$blend,
                    blend.threshold = opt$blend_threshold,
                    label = opt$label,
                    label.size = opt$label_size,
                    repel = opt$repel,
                    ncol = opt$ncol,
                    coord.fixed = opt$coord_fixed,
                    by.col = opt$do_not_by_col)
} else if ( opt$plot_type == 'VlnPlot' ) { 
plot_object <- VlnPlot(object = seurat_object,
                    features = features,
                    cols = cols,
                    pt.size = opt$pt_size,
                    idents = opt$idents,
                    sort = opt$sort,
                    assay = opt$assay,
                    group.by = group_by,
                    split.by = opt$split_by,
                    adjust = opt$adjust,
                    y.max = opt$y_max,
                    same.y.lims = opt$same_y_lims,
                    log = opt$log,
                    ncol = opt$ncol,
                    slot = opt$slot,
                    split.plot = opt$split_plot,
                    stack = opt$stack,
                    fill.by = opt$fill_by,
                    flip = opt$flip)
} else if ( opt$plot_type == 'DotPlot' ) { 
plot_object <- DotPlot(object = seurat_object,
                    assay = opt$assay,
                    features = features,
                    cols = cols_dot_plot,
                    col.min = opt$col_min,
                    col.max = opt$col_max,
                    dot.min = opt$dot_min,
                    dot.scale = opt$dot_scale,
                    idents = idents,
                    group.by = group_by,
                    split.by = opt$split_by,
                    cluster.idents = opt$cluster_idents,
                    scale = opt$do_not_scale,
                    scale.by = opt$scale_by,
                    scale.min = opt$scale_min,
                    scale.max = opt$scale_max)
}
if (!is.null(opt$output_rds_file)) {

saveRDS(file = opt$output_rds_file,
                    object = plot_object)
}

ggsave(filename = opt$plot_out,
                    plot = plot_object,
                    device = opt$plot_format,
                    scale = opt$scale_factor,
                    width = opt$width,
                    height = opt$height,
                    units = opt$units,
                    dpi = opt$dpi,
                    limitsize = opt$do_not_limitsize,
                    bg = opt$bg)
