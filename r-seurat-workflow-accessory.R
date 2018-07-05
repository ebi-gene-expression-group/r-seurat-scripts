# Wrap optparse's functionality to add support for mandatory arguments etc

rsw_parse_args <- function(option_list, mandatory=c()){
  
  parser = OptionParser(option_list=option_list)
  opt = parse_args(parser, convert_hyphens_to_underscores = TRUE)
  
  # Check options
  
  for (comp in mandatory){
    if (is.na(opt[[comp]])){
      print_help(parser)
      stop(paste0("Mandatory argment '", comp, "' not supplied"))
    }
  }
  
  opt
}

# Check metadata variables against Seurat object

check_metadata <- function(seurat_object, varnames){
  for ( vn in varnames ){
    if (! vn %in% colnames(seurat_object@meta.data) ){
      stop(paste(paste0("'", vn, "'"), 'not a valid metadata variable for this object'))
    }
  }
}

# Check specified cells against object

check_cells <- function(cell_names, seurat_object){
  if (! all(cell_names %in% seurat_object@cell.names)){
    stop("Some specified cells not present in Seurat object")
  }
}

# Parse numeric parameters and make sure vectors are the right length

parse_numeric <- function(opt, varname, val_for_na=NA, length=1){
  if (is.na(opt[[varname]])){
    return(rep(val_for_na, length))
  }else{
    vals <- split_string(opt[[varname]])
  }
  
  if ( any(c(grepl('[^0-9.]', vals)))){
    stop("Non-numeric filters supplied")
  }else{
    return(as.numeric(vals))
  }
}

# Split a string to a character vector

split_string <- function(x, sep=','){
  unlist(strsplit(x, sep))
}