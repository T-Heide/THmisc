#' Get sample annotations from sample barcodes.
#'
#' @param barcodes A character vector of sample barcodes.
#' @param extract Logical flag if barcodes should be extracted from 'barcodes' (Default: False).
#'
#' @return A data frame with sample annotations, rownames are set if there are no duplicated barcodes. 
#' @export
#'
#' @examples annotation_from_barcode("EPICC_C516_Z1_B1_D1") # a EPICC sample barcode
annotation_from_barcode = function(barcodes, extract=FALSE) {
  
  # quick check of input
  checkmate::assertCharacter(barcodes, null.ok = FALSE)
  checkmate::assertFlag(extract)
  
  # try annotation with EPICC project barcodes
  try({
    annot = annotation_from_barcode_epicc(barcodes, extract)
    return(annot)
  }, silent = TRUE)
  
  stop("Couldn't annotate barcodes.\n")
}
