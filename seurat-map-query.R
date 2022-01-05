#!/usr/bin/env Rscript
# This script has been automatically generated through
#
# YAML2RScript.py -i ../r-seurat-scripts/cli-generated/manually_crafted_YAML/seurat-map-query.yaml -o ../r-seurat-scripts/seurat-map-query.R
#
# to change this file edit the input YAML and re-run the above command

suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(SeuratDisk))

option_list <- list(
    make_option(
        c("-q", "--query-object-file"),
        action = "store",
        metavar = "Input file",
        type = "character",
        help = "Query file with Seurat object in either RDS-Seurat, Loom or SCE"
    ),
    make_option(
        c("--query-format"),
        action = "store",
        default = "seurat",
        metavar = "Input format",
        type = "character",
        help = "Either loom, seurat, anndata or singlecellexperiment for the input format to read."
    ),
    make_option(
        c("-a", "--anchors-object-file"),
        action = "store",
        metavar = "Seurat transfer anchors file",
        type = "character",
        help = "Input file with Seurat object with anchors in either RDS-Seurat, Loom or SCE"
    ),
    make_option(
        c("--anchors-format"),
        action = "store",
        default = "seurat",
        metavar = "Anchors format",
        type = "character",
        help = "Either loom, seurat, anndata or singlecellexperiment for the anchors format to read."
    ),
    make_option(
        c("-r", "--reference-object-file"),
        action = "store",
        default = NULL,
        metavar = "Seurat reference object file",
        type = "character",
        help = "Input file with Seurat object with reference (and UMAP) in either RDS-Seurat, Loom or SCE"
    ),
    make_option(
        c("--reference-format"),
        action = "store",
        default = "seurat",
        metavar = "Anchors format",
        type = "character",
        help = "Either loom, seurat, anndata or singlecellexperiment for the reference format to read."
    ),
    make_option(
        c("--refdata-field-or-assay"),
        action = "store",
        type = "character",
        help = "The name of the metadata field or assay from the reference object provided. This requires the reference parameter to be specified."
    ),
    make_option(
        c("--new-reduction-name"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Name for new integrated dimensional reduction."
    ),
    make_option(
        c("--reference-reduction"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Name of reduction to use from the reference for neighbor finding"
    ),
    make_option(
        c("--reduction-model"),
        action = "store",
        default = NULL,
        type = "character",
        help = "DimReduc object name that contains the umap model"
    ),
    make_option(
        c("--transferdata-args"),
        action = "store",
        default = NULL,
        type = "character",
        help = "A named list of additional arguments to TransferData, written in R syntax .ie list( argument = 'value' )"
    ),
    make_option(
        c("--integrated-embedding-args"),
        action = "store",
        default = NULL,
        type = "character",
        help = "A named list of additional arguments to IntegrateEmbeddings, written in R syntax .ie list( argument = 'value' )"
    ),
    make_option(
        c("--project-umap-args"),
        action = "store",
        default = NULL,
        type = "character",
        help = "A named list of additional arguments to ProjectUMAP, written in R syntax .ie list( argument = 'value' )"
    ),
    make_option(
        c("-o", "--output-object-file"),
        action = "store",
        type = "character",
        help = "FILE IN"
    ),
    make_option(
        c("--output-format"),
        action = "store",
        default = "seurat",
        type = "character",
        help = "Either loom, seurat, anndata or singlecellexperiment for the anchors format to read.'"
    )
)

opt <- wsc_parse_args(option_list, 
                      mandatory = c("query_object_file", "anchors_object_file", "refdata_field_or_assay", "output_object_file"))
                

if (!file.exists(opt$query_object_file)) {
    stop((paste("File", opt$query_object_file, "does not exist")))
}



if (!file.exists(opt$anchors_object_file)) {
    stop((paste("File", opt$anchors_object_file, "does not exist")))
}



if (!file.exists(opt$reference_object_file)) {
    stop((paste("File", opt$reference_object_file, "does not exist")))
}



if (!is.null(opt$refdata_field_or_assay)) {
    opt$refdata_field_or_assay <- eval(parse(text = opt$refdata_field_or_assay))
}



if (!is.null(opt$transferdata_args)) {
    opt$transferdata_args <- eval(parse(text = opt$transferdata_args))
}



if (!is.null(opt$integrated_embedding_args)) {
    opt$integrated_embedding_args <- eval(parse(text = opt$integrated_embedding_args))
}



if (!is.null(opt$project_umap_args)) {
    opt$project_umap_args <- eval(parse(text = opt$project_umap_args))
}


load_seurat4_packages_for_format(formats = c(opt$query_format, opt$anchors_format, opt$reference_format))

seurat_object <- read_seurat4_object(input_path = opt$query_object_file,
                    format = opt$query_format)

anchors_object <- read_seurat4_object(input_path = opt$anchors_object_file,
                    format = opt$anchors_format)

reference_object <- read_seurat4_object(input_path = opt$reference_object_file,
                    format = opt$reference_format)

mapped_object <- MapQuery(query = seurat_object,
                    anchorset = anchors_object,
                    reference = reference_object,
                    refdata = opt$refdata_field_or_assay,
                    new.reduction.name = opt$new_reduction_name,
                    reference.reduction = opt$reference_reduction,
                    reduction.model = opt$reduction_model,
                    transferdata.args = opt$transferdata_args,
                    integrateembeddings.args = opt$integrated_embedding_args,
                    projectumap.args = opt$project_umap_args)

write_seurat4_object(seurat_object = mapped_object,
                    output_path = opt$output_object_file,
                    format = opt$output_format)
