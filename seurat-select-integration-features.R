#!/usr/bin/env Rscript
# This script has been automatically generated through
#
# YAML2RScript.py -i ../r-seurat-scripts/cli-generated/manually_crafted_YAML/seurat-select-integration-features.yaml -o ../r-seurat-scripts/seurat-select-integration-features.R
#
# to change this file edit the input YAML and re-run the above command

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(workflowscriptscommon))

option_list <- list(
    make_option(
        c("-i", "--input-object-files"),
        action = "store",
        metavar = "Input files",
        type = "character",
        help = "A list of Seurat objects between which to find anchors for downstream integration. Input files in either RDS-Seurat, Loom or SCE, comma separated."
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
        c("--nfeatures"),
        action = "store",
        default = 2000,
        metavar = "Number of features",
        type = "integer",
        help = "Number of features to return"
    ),
    make_option(
        c("--assay-list"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Name or vector of assay names (one for each object) from which to pull the variable features."
    ),
    make_option(
        c("--do-not-verbose"),
        action = "store_false",
        default = TRUE,
        type = "logical",
        help = "Print messages"
    ),
    make_option(
        c("--fvf-nfeatures"),
        action = "store",
        default = 2000,
        type = "integer",
        help = "nfeatures for FindVariableFeatures. Used if VariableFeatures have not been set for any object in input."
    ),
    make_option(
        c("--file-out"),
        action = "store",
        metavar = "Rdata file with features",
        type = "character",
        help = "FILE IN"
    )
)

opt <- wsc_parse_args(option_list, 
                      mandatory = c("input_object_files", "file_out"))
                
load_seurat4_packages_for_format(formats = c(opt$input_format))

seurat_objects <- read_multiple_seurat4_objects(input_path_list = opt$input_object_files,
                    format = opt$input_format)

features <- SelectIntegrationFeatures(object.list = seurat_objects,
                    nfeatures = opt$nfeatures,
                    assay = opt$assay_list,
                    verbose = opt$do_not_verbose,
                    fvf.nfeatures = opt$fvf_nfeatures)

saveRDS(object = features,
                    file = opt$file_out)
