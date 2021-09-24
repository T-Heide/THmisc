#' Get sample annotations from sample barcodes.
#'
#' @param barcodes A character (or factor) vector of sample barcodes.
#' @param extract Logical flag if barcodes should be extracted from 'barcodes' (Default: False).
#'
#' @return A data frame with sample annotations, rownames are set if there are no duplicated barcodes. 
#' @export
#'
#' @examples annotation_from_barcode("EPICC_C516_Z1_B1_D1") # a EPICC sample barcode
annotation_from_barcode = function(barcodes, extract=FALSE) {
  
  # quick check of input
  if (is.factor(barcodes)) barcodes = as.character(barcodes)
  checkmate::assertCharacter(barcodes, null.ok = FALSE)
  checkmate::assertFlag(extract)
  
  if (length(barcodes) == 0)
    return(NULL)
  
  # if any barcodes are duplicated only annotate unique barcodes
  if (any(duplicated(barcodes))) {
    annot = annotation_from_barcode(unique(barcodes))[barcodes,]
    rownames(annot) = NULL
    return(annot)
  }
  
  # try annotation with all known barcode fucntions
  annot_fun_to_test = 
    list(
      THmisc::annotation_from_barcode_epicc,
      THmisc::annotation_from_barcode_tcga
    )
  
  for (annot_fun in annot_fun_to_test) {
    try({
      annot = annotation_from_barcode_epicc(barcodes, extract)
      rownames(annot) = barcodes
      return(annot)
    }, silent = TRUE)
  }
  
  stop("Couldn't annotate barcodes.\n")
}
