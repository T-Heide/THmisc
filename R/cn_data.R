generate_segment_input = function(d) {

  if (inherits(d, "CollapsedVCF")) {
    sample_ids = colnames(d)
    chrs = as.character(GenomeInfoDb::seqnames(d))
    pos = BiocGenerics::start(d)
    mut_gr = GenomicRanges::GRanges(chrs, IRanges::IRanges(pos, pos))
    mut_ids = rownames(d)
  } else if (inherits(d, "GRanges")) {
    mut_gr = d
    mut_ids = as.character(mut_gr)
  } else if (is.character(d)) {
    mut_gr = as(gsub("_.*", "", d), "GRanges")
    mut_ids = d
  } else if (is.data.frame(d)) {
    if (all(c("chr", "start", "end") %in% colnames(d))){
      mut_gr = with(d, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
      mut_ids = as.character(mut_gr)
    } else {
      stop()
    }
  } else {
    stop()
  }


  return(list(mut_gr, mut_ids))
}

#' Get CN annotations from a list of cn data.
#'
#' @param mutation_data A CollapsedVCF object, GRanges object or a mutation id.
#' @param cna_segments A (named) list of CN segment data.
#' @param sample_ids Indices of samples from 'cna_segments' to include in the output.
#'
#' @return A matrix containing annotations for each mutation object.
#' @export
#'
get_cnas = function(mutation_data, cna_segments, sample_ids=NULL) {

  mdata =
    tryCatch({
      generate_segment_input(mutation_data)
      }, error = function(e) stop("Invalid 'mutation_data' input.\n")
    )

  mut_gr = mdata[[1]]
  mut_ids = mdata[[2]]


  if (is.null(sample_ids)) {
    stop("Missing 'sample_ids' please pass these as input.\n")
  }

  # check for missing segment data:
  sample_with_cnas = sample_ids %in% names(cna_segments)
  if (any(!sample_with_cnas)) {
    missing = paste0(sample_ids[!sample_with_cnas],collapse=", ")
    warning("Missing segment data: ", missing, ".\n", sep="")
  }

  # fill cna matrix:
  cna_status = matrix(NA, NROW(mut_gr), length(sample_ids))
  dimnames(cna_status) = list(mut_ids, sample_ids)

  for (i in which(sample_with_cnas)) { # each sample

    # get samples cna data convert to granges
    sample_cna = cna_segments[[sample_ids[i]]] # samples segment data
    wh_cols = c(chr="chromosome",start="start.pos",end="end.pos")

    sample_cna_gr =
      magrittr::set_colnames(sample_cna[, wh_cols], names(wh_cols)) %>%
      as.data.frame() %>%
      as("GRanges")

    # find overlaps add cn to cna_status matrix
    ol = findOverlaps(mut_gr, sample_cna_gr, type="within", select="first")
    cna_status[,i] = c(sample_cna$CNt[ol])

    # try to calculate average of segments for missing data
    wh_missing = which(is.na(ol))
    if (length(wh_missing) == 0) next()

    ol = findOverlaps(mut_gr[wh_missing,], sample_cna_gr, type="any")
    if (length(ol) == 0) next()

    cna_status[wh_missing,i] =
      vapply(wh_missing, function(i) {
        wh = queryHits(ol) == i
        if (sum(wh) == 0) return(as.numeric(NA))
        sw = pintersect(mut_gr[wh,], sample_cna_gr[subjectHits(ol[wh])])
        weighted.mean(sample_cna$CNt[[subjectHits(ol[wh])]], width(sw))
      }, numeric(1))
  }

  return(cna_status)
}

#' Get A&B allele numbers from a list of CN data.
#'
#' @param mutation_data A CollapsedVCF object, GRanges object or a mutation id.
#' @param cna_segments A (named) list of CN segment data that contain A and B allele data.
#' @param sample_ids Indices of samples from 'cna_segments' to include in the output.
#'
#' @return A matrix containing annotations for each mutation object.
#' @export
#'
get_ab_alleles = function(mutation_data, cna_segments, sample_ids=NULL) {

  if ("CollapsedVCF" %in% class(mutation_data) & is.null(sample_ids)) {
    sample_ids = colnames(mutation_data)
  }

  mutation_data = generate_segment_input(mutation_data)[[1]]
  chrs = as.character(GenomicRanges::seqnames(mutation_data))
  pos = GenomicRanges::start(mutation_data)

  if (!is.null(sample_ids)) {
    chrs = as.character(GenomicRanges::seqnames(mutation_data))
    pos = GenomicRanges::start(mutation_data)
  }


  # check for missing segment data:
  sample_with_cnas = sample_ids %in% names(cna_segments)
  if (any(!sample_with_cnas)) {
    missing = paste0(sample_ids[!sample_with_cnas], collapse=", ")
    warning("Missing segment data: ", missing, ".\n", sep="")
  }

  # fill cna matrix:
  cna_status = matrix(NA, NROW(mutation_data), length(sample_ids))
  dimnames(cna_status) = list(rownames(mutation_data), sample_ids)

  for (i in which(sample_with_cnas)) { # each sample
    sample_cna = cna_segments[[sample_ids[i]]] # samples segment data

    for (j in seq_len(NROW(sample_cna))) { # each cna segment
      wh = chrs == sample_cna$chromosome[j] &
        between(pos, sample_cna$start.pos[j], sample_cna$end.pos[j])

      if (is.na(sample_cna$A[j]) | is.na(sample_cna$B[j])) {
        ab_status = NA
      } else {
        a_count = rep("A", sample_cna$A[j])
        b_count = rep("B", sample_cna$B[j])
        ab_status = paste0(c(a_count, b_count), collapse="")
      }
      cna_status[wh,i] = ab_status
    }
  }

  return(cna_status)
}

