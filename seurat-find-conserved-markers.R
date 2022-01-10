#!/usr/bin/env Rscript
# This script has been automatically generated through
#
# YAML2RScript.py -i ../r-seurat-scripts/cli-generated/manually_crafted_YAML/seurat-find-conserved-markers.yaml -o ../r-seurat-scripts/seurat-find-conserved-markers.R
#
# to change this file edit the input YAML and re-run the above command

suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(SeuratDisk))
suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(Seurat))

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
        c("--ident-1"),
        action = "store",
        default = "",
        type = "character",
        help = "Identity class to define markers for"
    ),
    make_option(
        c("--ident-2"),
        action = "store",
        default = NULL,
        type = "character",
        help = "A second identity class for comparison. If NULL (default) - use all other cells for comparison."
    ),
    make_option(
        c("--grouping-var"),
        action = "store",
        default = "",
        type = "character",
        help = "grouping variable"
    ),
    make_option(
        c("--assay"),
        action = "store",
        default = "RNA",
        type = "character",
        help = "of assay to fetch data for (default is RNA)"
    ),
    make_option(
        c("--slot"),
        action = "store",
        default = "data",
        type = "character",
        help = "Slot to pull data from; note that if test.use is negbinom, poisson, or DESeq2, slot  will be set to counts."
    ),
    make_option(
        c("--meta-method"),
        action = "store",
        default = "metap::minimump",
        type = "character",
        help = "method for combining p-values. Should be a function from the metap package (NOTE: pass the function, not a string)"
    ),
    make_option(
        c("--reduction"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Reduction to use in differential expression testing - will test for DE on cell embeddings"
    ),
    make_option(
        c("--features"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Genes to test. Default is to use all genes"
    ),
    make_option(
        c("--logfc-threshold"),
        action = "store",
        default = 0,
        type = "integer",
        help = "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing logfc.threshold speeds up the function, but can miss weaker signals."
    ),
    make_option(
        c("--test-use"),
        action = "store",
        default = "wilcox",
        type = "character",
        help = "Identifies differentially expressed genes between two groups using (see options)"
    ),
    make_option(
        c("--min-pct"),
        action = "store",
        default = 0,
        type = "integer",
        help = "only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1"
    ),
    make_option(
        c("--min-diff-pct"),
        action = "store",
        default = "-Inf",
        type = "double",
        help = "only test genes that show a minimum difference in the fraction of detection between the two groups. Set to -Inf by default"
    ),
    make_option(
        c("--do-not-verbose"),
        action = "store_false",
        default = TRUE,
        type = "logical",
        help = "Print a progress bar once expression testing begins"
    ),
    make_option(
        c("--only-pos"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Only return positive markers (FALSE by default)"
    ),
    make_option(
        c("--max-cells-per-ident"),
        action = "store",
        default = "Inf",
        type = "double",
        help = "Down sample each identity class to a max number. Default is no downsampling. Not activated by default (set to Inf)"
    ),
    make_option(
        c("--random-seed"),
        action = "store",
        default = 1,
        type = "integer",
        help = "Random seed for downsampling"
    ),
    make_option(
        c("--latent-vars"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Variables to test, used only when test.use is one of LR, negbinom, poisson, or MAST."
    ),
    make_option(
        c("--min-cells-feature"),
        action = "store",
        default = 3,
        type = "integer",
        help = "Minimum number of cells expressing the feature in at least one of the two groups, currently only used for poisson and negative binomial tests"
    ),
    make_option(
        c("--min-cells-group"),
        action = "store",
        default = 3,
        type = "integer",
        help = "Minimum number of cells in one of the groups"
    ),
    make_option(
        c("--pseudocount-use"),
        action = "store",
        default = 1,
        type = "integer",
        help = "Pseudocount to add to averaged expression values when calculating logFC. 1 by default."
    ),
    make_option(
        c("--mean-fxn"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Function to use for fold change or average difference calculation. If NULL, the appropriate function will be chose according to the slot used"
    ),
    make_option(
        c("--fc-name"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Name of the fold change, average difference, or custom function column in the output data.frame. If NULL, the fold change column will be named according to the logarithm base (eg, avg_log2FC), or if using the scale.data slot avg_diff."
    ),
    make_option(
        c("--base"),
        action = "store",
        default = 2,
        type = "integer",
        help = "The base with respect to which logarithms are computed."
    ),
    make_option(
        c("-o", "--markers_output_file"),
        action = "store",
        type = "character",
        help = "Output path for tab separated conserved marker genes file."
    )
)

opt <- wsc_parse_args(option_list,
                      mandatory = c("input_object_file", "markers_output_file"))


if (!file.exists(opt$input_object_file)) {
    stop((paste("File", opt$input_object_file, "does not exist")))
}


load_seurat4_packages_for_format(formats = c(opt$query_format, opt$anchors_format, opt$reference_format))

seurat_object <- read_seurat4_object(input_path = opt$input_object_file,
                    format = opt$input_format)

conserved_markers <- FindConservedMarkers(object = seurat_object,
                    ident.1 = opt$ident_1,
                    ident.2 = opt$ident_2,
                    grouping.var = opt$grouping_var,
                    assay = opt$assay,
                    slot = opt$slot,
                    meta.method = eval(parse(text=opt$meta_method)),
                    reduction = opt$reduction,
                    features = opt$features,
                    logfc.threshold = opt$logfc_threshold,
                    test.use = opt$test_use,
                    min.pct = opt$min_pct,
                    min.diff.pct = opt$min_diff_pct,
                    verbose = opt$do_not_verbose,
                    only.pos = opt$only_pos,
                    max.cells.per.ident = opt$max_cells_per_ident,
                    random.seed = opt$random_seed,
                    latent.vars = opt$latent_vars,
                    min.cells.feature = opt$min_cells_feature,
                    min.cells.group = opt$min_cells_group,
                    pseudocount.use = opt$pseudocount_use,
                    mean.fxn = opt$mean_fxn,
                    fc.name = opt$fc_name,
                    base = opt$base)

write.table(x = conserved_markers,
                    file = opt$markers_output_file,
                    sep = "\t")
