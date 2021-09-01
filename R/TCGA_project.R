#' TCGA project: Create a annotation for a vector of TCGA barcodes.
#'
#' @param barcodes A character vector of TCGA sample barcodes.
#' @param extract Logical flag if barcodes should be extracted from 'barcodes' (Default: False).
#'
#' @return A data frame with sample annotations.
#' @export
#'
#' @examples annotation_from_barcode_tcga("TCGA-02-0001-01C-01D")
annotation_from_barcode_tcga = function(barcodes, extract=TRUE) {

  # Wildcard identifying TCGA barcode specifications:
  # TCGA-TSS-Patient-Sample+Vial-Portion+Analyte-Plate-Center
  # TCGA-00-AAAA-00A-00A-0000-00
  wildcard = paste0(
    "TCGA", "-[[:alnum:]]{2}", "-[[:alnum:]]{4}",
    "(-[[:digit:]]{2}","([[:upper:]]{1}", "-([[:digit:]]{2}",
    "([DGHRTWX]{1}", "(-[[:alnum:]]{4}", "(-[[:digit:]]{2}",
    ")?)?)?)?)?)?"
  )

  # Stop if some barcodes do not meet the specs:
  err = sum(!grepl(wildcard, barcodes))
  if (err) {
    stop(sprintf("Error: %d barcodes do not meet EPICC specification.\n", err))
  }

  # Extract barcodes from strings (e.g. file names):
  if (extract) {
    barcodes_extracted = regmatches(barcodes, regexpr(wildcard, barcodes))
    n_fixed = sum(barcodes_extracted != barcodes)
    if (n_fixed) {
      if (extract) {
        msg = sprintf("Warning: Extracted %d barcodes from supplied strings.\n", n_fixed)
        warning(msg)
      } else {
        stop(sprintf("Error: %d barcodes do not meet EPICC specification.\n", n_fixed))
      }
    }
  }

  # Create annotation data.frame from barcodes:
  splt_name = strsplit(barcodes_extracted, "-")
  splt_data = data.frame(do.call("rbind", lapply(splt_name, "[", 1:7)))
  cnames = c("Study","TSS.Code","Patient","C4","C5","Plate","Center")
  colnames(splt_data) = cnames
  rownames(splt_data) = barcodes
  splt_data$Barcode = barcodes

  # Split some elements:
  splt_data$Type.ID = substr(splt_data$C4, 1, 2)
  splt_data$Vial    = substr(splt_data$C4, 3, 3)
  splt_data$Portion = substr(splt_data$C5, 1, 2)
  splt_data$Analyte = substr(splt_data$C5, 3, 3)
  splt_data[splt_data == ""] <- NA

  # Get additional annotations from code tables:
  splt_data = cbind(splt_data, tissueSourceSite[splt_data$TSS.Code,-1])
  splt_data = cbind(splt_data, centerCode[splt_data$Center,-1])
  splt_data = merge(splt_data, diseaseStudy, by="Study.Name", all.x=TRUE)
  splt_data = cbind(splt_data, portionAnalyte[splt_data$Analyte, -1, drop=0])
  id_type   = as.character(as.numeric(splt_data$Type.ID))
  splt_data = cbind(splt_data, sampleType[id_type, -1, drop=0])

  # Isolate selected columns and update the column names:
  sel_cols = c("Barcode"="Barcode",
                "Project.ID"="Study",
                "TSS.Code"="TSS.Code",
                "Study.Abbreviation"="Study.Abbreviation",
                "Study.Name"="Study.Name",
                "Patient.ID"="Patient",
                "Type.ID"="Type.ID",
                "Type.Short.Letter.Code"="Short.Letter.Code",
                "Type.Definition"="Definition",
                "Vial"="Vial",
                "Portion"="Portion",
                "Analyte.ID"="Analyte",
                "Analyte.Name"="Definition",
                "Plate.ID"="Plate",
                "Center.ID"="Center",
                "Center.Name"="Center.Name",
                "Center.Type"="Center.Type",
                "Center.Display.Name"="Display.Name",
                "Source.ID"="BCR",
                "Source.Site"="Source.Site")

  splt_data_ret = splt_data[,sel_cols]
  colnames(splt_data_ret) = names(sel_cols)
  rownames(splt_data_ret) = NULL
  
  return(splt_data_ret)
}
