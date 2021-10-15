#' Convert long peptide names to short ones.
#'
#' @param annot A character vector of annotations with peptide abbreviations
#'
#' @return Also a character vector, put with short peptide names
#' @export
#'
#' @examples long_to_short_peptide_names(c("TerAlaVal"))
long_to_short_peptide_names =  function(annot) {
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1867422/
  
  annot = as.character(annot)
  
  liftv = c("=", "*", seqinr::a())
  names(liftv) = c("%3D", "Ter", seqinr::aaa()) 
  
  for (i in seq_along(liftv)){
    annot = gsub(names(liftv)[i], liftv[i], annot, fixed=TRUE)
  }
  
  return(annot)
}

