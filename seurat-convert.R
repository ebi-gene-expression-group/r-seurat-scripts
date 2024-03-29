#!/usr/bin/env Rscript
# This script has been automatically generated through
#
# YAML2RScript.py -i ../r-seurat-scripts/cli-generated/manually_crafted_YAML/seurat-convert.yaml -o ../r-seurat-scripts/seurat-convert.R
#
# to change this file edit the input YAML and re-run the above command

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(SeuratDisk))

option_list <- list(
    make_option(
                    c("-i", "--input-object-file"),
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
                    help = "Either loom, seurat, h5seurat, anndata or singlecellexperiment for the input format to read.")
,
    make_option(
                    c("-o", "--output-object-file"),
                    action = "store",
                    type = "character",
                    help = "Path to the output file, when using Loom as output, the final file will contain .loom at the end of the given file name.")
,
    make_option(
                    c("--output-format"),
                    action = "store",
                    default = "seurat",
                    type = "character",
                    help = "Either seurat, h5seurat, loom, singlecellexperiment or h5seurat (partial support)")
,
    make_option(
                    c("--input-assay"),
                    action = "store",
                    default = "RNA",
                    type = "character",
                    help = "The assay to be read from the input object, if it corresponds")
)

opt <- wsc_parse_args(option_list,
                      mandatory = c("input_object_file", "output_object_file"))


if (!file.exists(opt$input_object_file)) {
    stop((paste("File", opt$input_object_file, "does not exist")))
}

if(opt$input_assay == "NULL") {
    opt$input_assay <- NULL
}

seurat_object <- read_seurat4_object(input_path = opt$input_object_file,
                    format = opt$input_format, assay = opt$input_assay)

write_seurat4_object(seurat_object = seurat_object,
                    output_path = opt$output_object_file,
                    format = opt$output_format)
