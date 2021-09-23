#!/usr/bin/env Rscript
# This script has been automatically generated through
#
# YAML2RScript.py -i ../r-seurat-scripts/cli-generated/manually_crafted_YAML/seurat-convert.yaml -o ../r-seurat-scripts/seurat-convert.R
#
# to change this file edit the input YAML and re-run the above command

suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(SeuratDisk))

option_list <- list(
    make_option(
                    c("-i", "--input-path"),
                    action = "store",
                    metavar = "Input file",
                    type = "character",
                    help = "Input file with Seurat object in either RDS-Seurat, Loom or SCE")
,
    make_option(
                    c("--input-format"),
                    action = "store",
                    default = "seurat",
                    metavar = "Input format",
                    type = "character",
                    help = "Either loom, seurat, anndata or singlecellexperiment for the input format to read.")
,
    make_option(
                    c("-o", "--output-path"),
                    action = "store",
                    type = "character",
                    help = "Path to the output file, when using Loom as output, the final file will contain .loom at the end of the given file name.")
,
    make_option(
                    c("--output-format"),
                    action = "store",
                    default = "seurat",
                    type = "character",
                    help = "Either seurat, loom, singlecellexperiment or h5seurat (partial support)")
)

opt <- wsc_parse_args(option_list, 
                      mandatory = c("input_path", "output_path"))


if (!file.exists(opt$input_path)) {
    stop((paste("File", opt$input_path, "does not exist")))
}


seurat_object <- read_seurat4_object(input_path = opt$input_path,
                    format = opt$input_format)

write_seurat4_object(seurat_object = seurat_object,
                    output_path = opt$output_path,
                    format = opt$output_format)
