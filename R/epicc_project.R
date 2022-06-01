#' EPICC Project: Get sample annotations from sample barcodes.
#'
#' @param barcodes A character vector of EPICC sample barcodes.
#' @param extract Logical flag if barcodes should be extracted from 'barcodes' (Default: False).
#'
#' @return A data frame with sample annotations, rownames are set if there are no duplicated barcodes.
#' @export
#' @import dplyr
#'
#' @examples annotation_from_barcode_epicc("EPICC_C516_Z1_B1_D1")
#' @examples annotation_from_barcode_epicc("Text EPICC_C516_Z1_B1_D1 Text", extract = TRUE)
annotation_from_barcode_epicc =
  function (barcodes, extract = FALSE)
  {

    epicc_regex =
      paste0(
        "(?P<read_barcode>",
        "(?P<lane_barcode>",
        "(?P<iteration_barcode>",
        "(?P<sample_barcode>",
        "(?P<project>EPICC)_",
        "(?P<patient>C[0-9]+)_",
        "(?P<region>[A-Z])(?P<region_number>[0-9]+)_",
        "(?P<sample_type>[BGL])(?P<sample_number>[0-9]+)_",
        "(?P<analyte>[DRCBL])(?P<analyte_number>[0-9]+))",
        "(?:_I(?P<iteration>[0-9]+))?)",
        "(?(iteration)_(?P<lane_id>L[0-9A-Za-z]+)|)?)",
        "(?(lane_id)_R(?P<read_number>[012]))?)"
      )

    numeric_columns =    # list of columns/elements that
      c("region_number", # should be converted to numeric values
        "sample_number",
        "analyte_number",
        "iteration",
        "read_number")


    # mapping of ids to annotations:
    analyte_id_names = c(
      "D" = "WGS",
      "L" = "LP-WGS",
      "C" = "ATAC-seq",
      "R" = "RNA-seq",
      "B" = "Bisulfit WGS"
    )

    type_id_names = c(
      "B" = "bulk",
      "G" = "gland",
      "L" = "interglandular",
      "Z" = "blood"
    )

    tt_order = c(
      "normal",
      "adenoma",
      "cancer"
    )

    if (!exists("msi_positiv")) { # default value for msi positive cases if variable does not exist globally
      msi_positiv = c("C536","C548","C516","C518","C552","C562")
    }

    if (!exists("msi_positiv_adenoma")) { # default value for msi positive cases if variable does not exist globally
      msi_positiv_adenoma = c("C516")
    }

    # check if the input is valid:
    if (!is.vector(barcodes)) {
      stop(paste("Argument", sQuote("barcodes"), "has to be a vector.\n"))
    } else {
      barcodes = as.character(barcodes)
    }

    if (!is.logical(extract)) {
      stop(paste("Argument", sQuote("extract"), "has to be a boolean.\n"))
    }

    # check for non matching barcodes:
    regexpr_result = regexpr(epicc_regex, barcodes, perl = TRUE)
    nerr = sum(attr(regexpr_result, "match.length") == -1, na.rm=TRUE)
    if (nerr) {
      stop(sprintf("Error: %d barcode(s) do not meet EPICC specification.\n", nerr))
    }


    # check if a valid barcode can be extracted from the input:
    barcodes_extracted = regmatches(barcodes, regexpr_result)
    n_extr = sum(barcodes_extracted != barcodes)
    if (n_extr) {
      if (extract) {
        msg = sprintf("Extracted %d barcode(s) from supplied strings.\n", n_extr)
        warning(msg)
      }
      else {
        msg = sprintf("Error: %d barcode(s) do not meet EPICC specification.\n", n_extr)
        stop(msg)
      }
      regexpr_result = regexpr(epicc_regex, barcodes_extracted,  perl=TRUE)
    }


    # get the annotation elements:
    annotation =
      regcapturedmatches(barcodes_extracted, regexpr_result) %>%
      data.frame(stringsAsFactors = FALSE) %>%
      dplyr::mutate(lane_barcode=ifelse(lane_id == "", NA, lane_barcode)) %>%
      dplyr::mutate(iteration_barcode=ifelse(iteration == "", NA, iteration_barcode)) %>%
      dplyr::mutate(read_barcode=ifelse(read_number == "", NA, read_barcode)) %>%
      dplyr::mutate(tissue_barcode=gsub("_[DRCBL][0-9]+$", "", sample_barcode))

    if (sum(duplicated(barcodes)) == 0) {
      rownames(annotation) = barcodes
    }


    # insert tissue type:
    annotation$tissue_type =
      with(annotation, {
        dplyr::case_when( # some exceptions from the rule ...
          region %in% c("F") & patient == "C542" ~ "cancer",
          region %in% c("C", "D") & patient == "C516" ~ "adenoma",
          region %in% c("E", "Z", "W", "S") ~ "normal",
          region %in% c("A", "B", "C", "D") ~ "cancer",
          region %in% c("F", "G", "H", "I") ~ "adenoma"
        ) %>% factor(tt_order, ordered=TRUE)
      })


    # insert long name for analytes:
    annotation$analyte_name =
      with(annotation, {
        analyte_id_names[as.character(analyte)] %>%
          factor(analyte_id_names, ordered=TRUE)
      })


    # insert long name for sample type:
    annotation$sample_type_name =
      with(annotation, {
        dplyr::case_when(
          region == "Z" ~ "blood",
          region == "S" ~ "blood",
          TRUE ~ type_id_names[as.character(sample_type)]
        ) %>% factor(c(type_id_names), ordered=TRUE)
      })


    # convert some cols to numeric:
    for (col in numeric_columns) {
      annotation[, col] = as.numeric(as.character(annotation[, col]))
    }


    # insert a label for each tumour (e.g. independed adenomas):
    group = paste0(annotation$patient, ".", annotation$tissue_type)
    wh_adenoma = grepl("adenoma", group) # add adenoma number to labels
    adenoma_regions = annotation$region[wh_adenoma]
    adenoma_regions = gsub("[CD]", "C+D", adenoma_regions)

    adenoma_region_label_list =
      split(adenoma_regions, annotation$patient[wh_adenoma]) %>%
      lapply(function(x) {
        xu = unique(x)
        if (length(xu) > 1) { l = paste0(" (", xu, ")") } else { l = "" }
        names(l) = xu
        return(l)
      }) %>% unlist()

    key_label = paste0(annotation$patient[wh_adenoma], ".", adenoma_regions)
    adenoma_labels = adenoma_region_label_list[key_label]
    group[wh_adenoma] = paste0(group[wh_adenoma], adenoma_labels)
    annotation$tumour_id = group


    # add msi status
    annotation$msi_status =
      with(annotation, {
        dplyr::case_when(
          tissue_type == "normal" ~ as.character(NA),
          tissue_type == "cancer" & patient %in% msi_positiv ~ "MSI",
          tissue_type == "adenoma" & patient %in% msi_positiv_adenoma ~ "MSI",
          TRUE ~ "MSS"
        )
      })

    return(annotation)
  }
