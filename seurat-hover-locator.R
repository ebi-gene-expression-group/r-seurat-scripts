#!/usr/bin/env Rscript
# This script has been automatically generated through
#
# YAML2RScript.py -i ../r-seurat-scripts/cli-generated/manually_crafted_YAML/seurat-hover-locator.yaml -o ../r-seurat-scripts/seurat-hover-locator.R
#
# to change this file edit the input YAML and re-run the above command

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(plotly))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(htmlwidgets))
suppressPackageStartupMessages(require(optparse))

option_list <- list(
    make_option(
        c("--plot-rds"),
        action = "store",
        metavar = "RDS with ggplot2 object from Seurat 4",
        type = "character",
        help = "FILE IN"
    ),
    make_option(
        c("--information-table"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Table with information for overlaying on the plot. Usually the result of calling FetchData on the original Seurat object."
    ),
    make_option(
        c("--do-not-stringsAsFactors"),
        action = "store_false",
        default = TRUE,
        type = "logical",
        help = "Load strings as factors."
    ),
    make_option(
        c("--key"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Key (index) for the information table"
    ),
    make_option(
        c("--do-not-keepLeadingZeros"),
        action = "store_false",
        default = TRUE,
        type = "logical",
        help = "If TRUE a column containing numeric data with leading zeros will be read as character, otherwise leading zeros will be removed and converted to numeric."
    ),
    make_option(
        c("--do-not-axes"),
        action = "store_false",
        default = TRUE,
        type = "logical",
        help = "Display x- and y-axes"
    ),
    make_option(
        c("--dark-theme"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Plot using a dark theme?"
    ),
    make_option(
        c("--output-html"),
        action = "store",
        type = "character",
        help = "HTML output file"
    ),
    make_option(
        c("--do-not-selfcontained"),
        action = "store_false",
        default = TRUE,
        type = "logical",
        help = "Whether to save the HTML as a single self-contained file (with external resources base64 encoded) or a file with external resources placed in an adjacent directory."
    ),
    make_option(
        c("--libdir"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Directory to copy HTML dependencies into (defaults to filename_files)."
    )
)

opt <- wsc_parse_args(option_list, 
                      mandatory = c("plot_rds", "output_html"))
                

if ( !file.exists(opt$plot_rds)) {
    stop((paste("File", opt$plot_rds, "does not exist")))
}



if (!is.null(opt$information_table) && !file.exists(opt$information_table)) {
    stop((paste("File", opt$information_table, "does not exist")))
}


ggplot2_object <- readRDS(file = opt$plot_rds)

if (!is.null(opt$information_table)) {

info_df <- fread(file = opt$information_table,
                    stringsAsFactors = opt$do_not_stringsAsFactors,
                    key = opt$key,
                    keepLeadingZeros = opt$do_not_keepLeadingZeros)
} else {
  info_df <- NULL
}

plotly_obj <- HoverLocator(plot = ggplot2_object,
                    information = info_df,
                    axes = opt$do_not_axes,
                    dark.theme = opt$dark_theme)

plotly_widget <- as.widget(x = plotly_obj)

saveWidget(widget = plotly_widget,
                    file = opt$output_html,
                    selfcontained = opt$do_not_selfcontained,
                    libdir = opt$libdir)
