# Wrap optparse's functionality to add support for mandatory arguments etc

rsw_parse_args <- function(option_list, mandatory=c()){
  
  parser = OptionParser(option_list=option_list)
  opt = parse_args(parser, convert_hyphens_to_underscores = TRUE)
  
  # Check options
  
  for (comp in mandatory){
    if (is.na(opt[[comp]])){
      print_help(parser)
      stop("One or more mandatory argments not supplied")
    }
  }
  
  opt
}