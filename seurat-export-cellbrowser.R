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
    c("--input-format"),
    action = "store",
    default = "seurat",
    type = 'character',
    help = "Either loom, seurat, anndata or singlecellexperiment for the input format to read."
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

suppressPackageStartupMessages(require(Seurat))
if(opt$input_format == "loom" ) {
  suppressPackageStartupMessages(require(SeuratDisk))
} else if(opt$input_format == "singlecellexperiment" ) {
  suppressPackageStartupMessages(require(scater))
}

# Input from serialized R object
seurat_object <- read_seurat4_object(input_path = opt$input_object_file, format = opt$input_format)

skip_markers = TRUE
if(!is.null(opt$markers_file)) {
  skip_markers = FALSE
}


# From https://raw.githubusercontent.com/satijalab/seurat-wrappers/cellbrowser/R/cellbrowser.R
# here in the meantime while this is not merged there.
ExportToCellbrowser <- function(
  object,
  dir,
  dataset.name = Project(object = object),
  reductions = "tsne",
  markers.file = NULL,
  cluster.field = "Cluster",
  cb.dir = NULL,
  port = NULL,
  skip.expr.matrix = FALSE,
  skip.metadata = FALSE,
  skip.reductions = FALSE,
  ...
) {
  vars <- c(...)
  if (is.null(x = vars)) {
    vars <- c("nCount_RNA", "nFeature_RNA")
    if (length(x = levels(x = Idents(object = object))) > 1) {
      vars <- c(vars, cluster.field)
      names(x = vars) <- c("", "", "ident")
    }
  }
  names(x = vars) <- names(x = vars) %||% vars
  names(x = vars) <- sapply(
    X = 1:length(x = vars),
    FUN = function(i) {
      return(ifelse(
        test = nchar(x = names(x = vars)[i]) > 0,
        yes = names(x = vars[i]),
        no = vars[i]
      ))
    }
  )
  if (!is.null(x = port) && is.null(x = cb.dir)) {
    stop("cb.dir parameter is needed when port is set")
  }
  if (!dir.exists(paths = dir)) {
    dir.create(path = dir)
  }
  if (!dir.exists(paths = dir)) {
    stop("Output directory ", dir, " cannot be created or is a file")
  }
  if (dataset.name == "SeuratProject") {
    warning("Using default project name means that you may overwrite project with the same name in the cellbrowser html output folder")
  }
  order <- colnames(x = object)
  enum.fields <- c()
  # Export expression matrix:
  if (!skip.expr.matrix) {
    # Relatively memory inefficient - maybe better to convert to sparse-row and write in a loop, row-by-row?
    df <- as.data.frame(x = as.matrix(x = GetAssayData(object = object)))
    df <- data.frame(gene = rownames(x = object), df, check.names = FALSE)
    gzPath <- file.path(dir, "exprMatrix.tsv.gz")
    z <- gzfile(gzPath, "w")
    message("Writing expression matrix to ", gzPath)
    write.table(x = df, sep = "\t", file = z, quote = FALSE, row.names = FALSE)
    close(con = z)
  }
  # Export cell embeddings
  embeddings.conf <- c()
  for (reduction in reductions) {
    if (!skip.reductions) {
      df <- Embeddings(object = object, reduction = reduction)
      if (ncol(x = df) > 2) {
        warning(
          'Embedding ',
          reduction,
          ' has more than 2 coordinates, taking only the first 2'
        )
        df <- df[, 1:2]
      }
      colnames(x = df) <- c("x", "y")
      df <- data.frame(cellId = rownames(x = df), df)
      fname <- file.path(dir, paste0(reduction, '.coords.tsv'))
      message("Writing embeddings to ", fname)
      write.table(
        x = df[order, ],
        sep = "\t",
        file = fname,
        quote = FALSE,
        row.names = FALSE
      )
    }
    conf <- sprintf(
      '{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}',
      reduction
    )
    embeddings.conf <- c(embeddings.conf, conf)
  }
  # Export metadata
  df <- data.frame(row.names = rownames(x = object[[]]))
  df <- FetchData(object = object, vars = names(x = vars))
  colnames(x = df) <- vars
  enum.fields <- Filter(
    f = function(name) {!is.numeric(x = df[[name]])},
    x = vars
  )
  if (!skip.metadata) {
    fname <- file.path(dir, "meta.tsv")
    message("Writing meta data to ", fname)
    df <- data.frame(Cell = rownames(x = df), df, check.names = FALSE)
    write.table(
      x = df[order, ],
      sep = "\t",
      file = fname,
      quote = FALSE,
      row.names = FALSE
    )
  }
  # Export markers
  markers.string <- ''
  if (!is.null(x = markers.file)) {
    ext <- file_ext(x = markers.file)
    fname <- paste0("markers.", ext)
    file.copy(from = markers.file, to = file.path(dir, fname))
    markers.string <- sprintf(
      'markers = [{"file": "%s", "shortLabel": "Seurat Cluster Markers"}]',
      fname
    )
  }
  config <- c(
    'name="%s"',
    'shortLabel="%1$s"',
    'exprMatrix="exprMatrix.tsv.gz"',
    '#tags = ["10x", "smartseq2"]',
    'meta="meta.tsv"',
    '# possible values: "gencode-human", "gencode-mouse", "symbol" or "auto"',
    'geneIdType="auto"',
    'clusterField="%s"',
    'labelField="%2$s"',
    'enumFields=%s',
    '%s',
    'coords=%s'
  )
  config <- paste(config, collapse = '\n')
  enum.string <- paste0(
    "[",
    paste(paste0('"', enum.fields, '"'), collapse = ", "),
    "]"
  )
  coords.string <- paste0(
    "[",
    paste(embeddings.conf, collapse = ",\n"),
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
  fname <- file.path(dir, "cellbrowser.conf")
  if (file.exists(fname)) {
    message(
      "`cellbrowser.conf` already exists in target directory, refusing to ",
      "overwrite it"
    )
  } else {
    cat(config, file = fname)
  }
  message("Prepared cellbrowser directory ", dir)
  if (!is.null(x = cb.dir)) {
    if (!py_module_available(module = "cellbrowser")) {
      stop(
        "The Python package `cellbrowser` is required to prepare and run ",
        "Cellbrowser. Please install it ",
        "on the Unix command line with `sudo pip install cellbrowser` (if root) ",
        "or `pip install cellbrowser --user` (as a non-root user). ",
        "To adapt the Python that is used, you can either set the env. variable RETICULATE_PYTHON ",
        "or do `require(reticulate) and use one of these functions: use_python(), use_virtualenv(), use_condaenv(). ",
        "See https://rstudio.github.io/reticulate/articles/versions.html; ",
        "at the moment, R's reticulate is using this Python: ",
        import(module = 'sys')$executable,
        ". "
      )
    }
    if (!is.null(x = port)) {
      port <- as.integer(x = port)
    }
    message("Converting cellbrowser directory to html/json files")
    cb <- import(module = "cellbrowser")
    cb$cellbrowser$build(dir, cb.dir)
    if (!is.null(port)) {
      message("Starting http server")
      cb$cellbrowser$stop()
      cb$cellbrowser$serve(cb.dir, port)
      Sys.sleep(time = 0.4)
      browseURL(url = paste0("http://localhost:", port))
    }
  }
  return(invisible(x = NULL))
}

ExportToCellbrowser(object = seurat_object, 
                    dir = opt$output_directory, 
                    markers.file = opt$markers_file, 
                    dataset.name = opt$study_name)
