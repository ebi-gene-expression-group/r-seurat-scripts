#!/usr/bin/env Rscript 

# This script currently calls an embedded function which is part of Seurat 3. As the package moves with Seurat 3, we will
# remove that function and rely on the native one.

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
    c("-o", "--output-directory"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'Seurat'.'"
  ),
  make_option(
    c("-n", "--study-name"),
    action= "store",
    default = "Seurat_study",
    type="character",
    help="Study name to be displayed in CellBrowser"
  ),
  make_option(
    c("-m", "--markers-file"),
    default = NULL,
    type="character",
    help="Path to markers file computed by Seurat"
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_directory'))
opt$study_name = gsub(" ","_", opt$study_name) # study names cannot have spaces.

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

# Input from serialized R object

seurat_object <- readRDS(opt$input_object_file)



# the following is copied from Seurat3's utilities.R and generics.R
# it was modified to work with Seurat2. At some point we probably should merge the two code bases.
# the markers "---" are used to mark the part that is used by the command line tool cbImportSeurat2

# ---
#' Export Seurat object for UCSC cell browser
#'
#' @param object Seurat object
#' @param dir output directory path
#' @param dataset.name name of the dataset
#' @param meta.fields vector of metadata fields to export
#' @param meta.fields.names list that defines metadata field names
#'                          after the export. Should map metadata
#'                          column name to export name
#' @param embeddings vector of embedding names to export
#' @param markers.file path to file with marker genes
#' @param cluster.field name of the metadata field containing cell cluster
#' @param cb.dir in which dir to create UCSC cell browser content root
#' @param port \emph{experimental} on which port to run UCSC cell browser after export
#'
#' @export
#'
#' @importFrom reticulate py_module_available
#' @importFrom reticulate import
#'
#' @examples
#' \dontrun{
#' ExportToCellbrowser(object = pbmc_small, dataset.name = "PBMC", dir = "out")
#' }
#'
ExportToCellbrowser <- function(
  object,
  dir,
  dataset.name,
  meta.fields = NULL,
  meta.fields.names = NULL,
  embeddings = c("tsne", "pca", "umap"),
  markers.file = NULL,
  cluster.field = NULL,
  port = NULL,
  cb.dir = NULL,
  markers.n = 100,
  skip.expr.matrix = FALSE,
  skip.markers = FALSE,
  all.meta = FALSE
) {
  if (!require("Seurat",character.only = TRUE)) {
    stop("This script requires that Seurat (V2 or V3) is installed")
  }
  message("Seurat Version installed: ", packageVersion("Seurat"))
  message("Object was created with Seurat version ", object@version)
  
  if (substr(object@version, 1, 1)!='2') {
    stop("can only process Seurat2 objects, version of rds is ", object@version)
  }
  
  #idents <- FetchData(object,c("ident"))$ident;
  idents <- object@ident # Idents() in Seurat3
  
  if (is.null(cluster.field)) {
    cluster.field = "Cluster"
  }
  
  if (is.null(meta.fields)) {
    meta.fields <- colnames(object@meta.data)
    if (length(levels(idents)) > 1) {
      meta.fields <- c(meta.fields, ".ident")
    }
  }
  if (!is.null(port) && is.null(cb.dir)) {
    stop("cb.dir parameter is needed when port is set")
  }
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  if (!dir.exists(dir)) {
    stop("Output directory ", dir, " cannot be created or is a file")
  }
  
  order <- object@cell.names
  enum.fields <- c()
  
  # Export expression matrix
  if (!skip.expr.matrix) { 
    mat <- as.matrix(object@raw.data)
    #mat <- as.matrix(GetAssayData(object = object, slot="counts"))
    df <- as.data.frame(mat, check.names=FALSE)
    df <- data.frame(gene=rownames(object@data), df, check.names = FALSE)
    gzPath <- file.path(dir, "exprMatrix.tsv.gz")
    z <- gzfile(gzPath, "w")
    message("Writing expression matrix to ", gzPath)
    write.table(x = df, sep="\t", file=z, quote = FALSE, row.names = FALSE)
    close(con = z)
  }
  
  # Export cell embeddings
  embeddings.conf <- c()
  dr <- object@dr
  for (embedding in embeddings) {
    #df <- Embeddings(object = object, reduction = embedding)
    emb <- dr[[embedding]]
    if (is.null(emb)) {
      message("Embedding ",embedding," does not exist in Seurat object. Skipping. ")
      next
    }
    #df <- slot(emb, "cell.embeddings")
    df <-  emb@cell.embeddings
    if (ncol(df) > 2) {
      warning('Embedding ', embedding, ' has more than 2 coordinates, taking only the first 2')
      df <- df[, 1:2]
    }
    colnames(df) <- c("x", "y")
    df <- data.frame(cellId = rownames(df), df, check.names=FALSE)
    fname <- file.path(
      dir,
      sprintf("%s.coords.tsv", embedding)
    )
    message("Writing embeddings to ", fname)
    write.table(df[order, ], sep="\t", file=fname, quote = FALSE, row.names = FALSE)
    conf <- sprintf(
      '{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}',
      embedding
    )
    embeddings.conf <- c(embeddings.conf, conf)
  }
  
  # Export metadata
  df <- data.frame(row.names = object@cell.names, check.names=FALSE)
  for (field in meta.fields) {
    if (field == ".ident") {
      df$Cluster <- idents
      enum.fields <- c(enum.fields, "Cluster")
    } else {
      name <- meta.fields.names[[field]]
      if (is.null(name)) {
        name <- field
      }
      #df[[name]] <- object[[field]][, 1]
      #df[[name]] <- FetchData(object, field)[, 1]
      df[[name]] <- object@meta.data[[field]]
      if (!is.numeric(df[[name]])) {
        enum.fields <- c(enum.fields, name)
      }
    }
  }
  df <- data.frame(Cell=rownames(df), df, check.names=FALSE)
  
  fname <- file.path(dir, "meta.tsv")
  message("Writing meta data to ", fname)
  write.table(df[order, ], sep="\t", file=fname, quote = FALSE, row.names=FALSE)
  
  # Export markers
  markers.string <- ''
  if (is.null(markers.file)) {
    ext <- "tsv"
  } else {
    ext <- tools::file_ext(markers.file)
  }
  file <- paste0("markers.", ext)
  fname <- file.path(dir, file)
  
  if (!is.null(markers.file) && !skip.markers) {
    message("Copying ", markers.file, " to ", fname)
    file.copy(markers.file, fname)
  }
  if (is.null(markers.file) && skip.markers) {
    file <- NULL
  }
  if (is.null(markers.file) && !skip.markers) {
    if (length(levels(idents)) > 1) {
      message("Running FindAllMarkers(), using wilcox test, min logfc diff 0.35, and writing top ", markers.n, ", cluster markers to ", fname)
      markers <- FindAllMarkers(object, do.print=TRUE, print.bar=TRUE, test.use="wilcox", logfc.threshold = 0.25)
      markers.helper <- function(x) {
        partition <- markers[x,]
        ord <- order(partition$p_val_adj < 0.05, -partition$avg_logFC)
        res <- x[ord]
        naCount <- max(0, length(x) - markers.n)
        res <- c(res[1:markers.n], rep(NA, naCount))
        return(res)
      }
      markers.order <- ave(rownames(markers), markers$cluster, FUN=markers.helper)
      top.markers <- markers[markers.order[!is.na(markers.order)],]
      write.table(top.markers, fname, quote=FALSE, sep="\t", col.names=NA)
    } else {
      message("No clusters found in Seurat object and no external marker file provided, so no marker genes can be computed")
      file <- NULL
    }
  }
  
  if (!is.null(file)) {
    markers.string <- sprintf(
      'markers = [{"file": "%s", "shortLabel": "Seurat Cluster Markers"}]',
      file
    )
  }
  
  config <- '
# This is a bare-bones, auto-generated cellbrowser config file.
# Look at https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/cellbrowser.conf
# for a full file that shows all possible options
name="%s"
shortLabel="%1$s"
exprMatrix="exprMatrix.tsv.gz"
#tags = ["10x", "smartseq2"]
meta="meta.tsv"
# possible values: "gencode-human", "gencode-mouse", "symbol" or "auto"
geneIdType="auto"
# file with gene,description (one per line) with highlighted genes, called "Dataset Genes" in the user interface
# quickGenesFile=quickGenes.csv
clusterField="%s"
labelField="%2$s"
enumFields=%s
%s
coords=%s'
enum.string <- paste0(
    "[",
    paste(paste0('"', enum.fields, '"'), collapse = ", "),
    "]"
  )
  coords.string <- paste0(
    "[",
    paste(embeddings.conf, collapse = ", "),
    "]"
  )
  config <- sprintf(
    config,
    dataset.name,
    cluster.field,
    enum.string,
    markers.string,
    coords.string
  )
  
  confPath = file.path(dir, "cellbrowser.conf")
  message("Writing cellbrowser config to ", confPath)
  cat(config, file=confPath)
  
  message("Prepared cellbrowser directory ", dir)
  
  
}

skip_markers = TRUE
if(!is.null(opt$markers_file)) {
  skip_markers = FALSE
}

ExportToCellbrowser(object = seurat_object, dir = opt$output_directory, skip.markers = skip_markers, markers.file = opt$markers_file, dataset.name = opt$study_name)
