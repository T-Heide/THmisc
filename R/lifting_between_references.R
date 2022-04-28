#' Internal function used to lift mutation data between genome versions.
#'
#' @param d Input data to work on. Options are i) a list of mutation vectors, ii) a GRanges object.
#' @param chain A chain file for the lifting.
#' @param geno A reference genome file of class 'BSgenome'.
#' @param check_new_alleles Logical flag indicating if match of reference alleles should be checked.
#' @param multi_summary_function Function used to summarise multiple segments in the new genome.
#'
#' @return A object with a structure similar to d with updated position data. 
#' @import GenomicRanges
#'
#' @examples 
lift_data_wrapper = function(d, chain, geno, check_new_alleles=TRUE, multi_summary_function=NULL) {
  
  checkmate::assertClass(geno, "BSgenome", null.ok = !check_new_alleles)
  checkmate::assert_file_exists(chain)
  checkmate::assert_flag(check_new_alleles)
  checkmate::assert_function(multi_summary_function, null.ok = TRUE)
  
  # convert data 'd' to a acceptable input:
  if (is.factor(d)) {
    d = as.character(d)
  }
  
  if (is.character(d)) {
    if (all(grepl("^(chr)[0-9XY]+:[0-9]+_[ACGT]+/[ACGT]+$", d))) {
      # mutation ids
      d = d %>%
        strsplit("[:_/]") %>%
        do.call(what=rbind) %>%
        data.frame(stringsAsFactors=FALSE) %>%
        magrittr::set_colnames(c("chr","pos","ref","alt"))
      
    } else if (grepl("^(chr)[0-9XY]+:[0-9]+-[0-9]+$", d)) {
      # position intervals 
      d = as(d, "GRanges")
      
    } else {
      # something else ...
      err_msg = 
        paste0(
          "Invalid or mixed string input for 'd'.\n\n",
          "  Valid formats are: \n",
          "   - '^(chr)[0-9]+:[0-9]+_[ACGT]+/[ACGT]+$' (mutation id)\n",
          "   - '^(chr)[0-9]+:[0-9]+-[0-9]+$' (position interval)"
        )
      
      stop(err_msg)
    }
  }
  
  # convert data to a granges object:
  # 1) A GRanges object.
  # 2) A mutation data.frame containing chr, pos, ref, alt columns
  if (is(d, "GRanges"))  {
    cur = d
  } else if (is.data.frame(d)){
    stopifnot(c("chr","pos","ref","alt") %in% colnames(d))
    iranges = IRanges::IRanges(start=as.numeric(d$pos), width=nchar(d$ref))
    cur = GenomicRanges::GRanges(d$chr, iranges)
  } else {
    stop("Invalid input for data 'd'\n")
  }
  
  # lift granges object to new genome
  ch = rtracklayer::import.chain(chain)
  new = rtracklayer::liftOver(cur, ch)
  
  # summary of multiple hits in genome, then drop multi hits:
  if (!is.null(multi_summary_function)) {
    new = multi_summary_function(new)
  }
  
  wh = sapply(new, NROW) == 1
  n_multi_hits = sum(!wh) > 0
  
  if (n_multi_hits) {
    warning(sprintf("Removed %d variant(s) missing in new genome or with multiple hits during liftover.\n", n_multi_hits))
  }
  
  d = d[wh,]
  cur = cur[wh,]
  new = unlist(new[wh])
  
  
  #
  if (is.data.frame(d)) {
    # is a data.frame (with chr, pos, ref and alt columns)
    # update these
    
    # get new ref allele:
    ref = as.character(Biostrings::getSeq(geno, new))
    
    # check ref allele match:
    if (check_new_alleles) {
      wh = ref == d$ref
      
      n_ref_missmatch = sum(!wh) > 0
      if (n_ref_missmatch) {
        warning(sprintf("Removed %d variants with new ref allele during liftover.\n", n_ref_missmatch))
      }
      
      d = d[wh,]
      cur = cur[wh,]
      new = new[wh,]
      ref = ref[wh]
    }
    
    # replace in old data:
    d$pos_old = d$pos
    d$pos = start(new)
    d$ref_old = d$ref
    d$ref = as.character(ref)
    d$chr_old = d$chr
    d$chr = as.character(GenomeInfoDb::seqnames(new))
    
  } else {
    # granges, copy over metadata columns
    S4Vectors::mcols(new) = S4Vectors::mcols(d)
    d = new
  }
  
  return(d)
}  


#' Function to lift mutations or objects between genome versions
#'
#' @param d Input data to work on. Options are i) a list of mutation ids, ii) a GRanges object.
#' @param ... Additional arguments passed to lift_data_wrapper (see ?lift_data_wrapper for details)
#'
#' @return A object with a structure similar to d with updated position data. 
#' @export
#'
#' @examples lift_hg38_to_hg19("chr2:148654004_T/C")
lift_hg38_to_hg19 = function(d, ...) {
  
  geno = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  
  chain = system.file(
    "extdata",
    "hg38ToHg19.over.chain",
    package = "THmisc",
    mustWork = TRUE
  )
  
  chain = getOption("chain_hg38_to_hg19") 
  if (is.null(chain)) 
    stop(paste0("Please set the ", sQuote("chain_hg38_to_hg19"), " option.\n"))
  
  lift_data_wrapper(d, chain, geno, ...)
}


#' Function to lift mutations or objects between genome versions
#'
#' @param d Input data to work on. Options are i) a list of mutation ids, ii) a GRanges object.
#' @param ... Additional arguments passed to lift_data_wrapper (see ?lift_data_wrapper for details)
#'
#' @return A object with a structure similar to d with updated position data. 
#' @export
#'
#' @examples lift_hg19_to_hg38("chr2:148654004_T/C")
lift_hg19_to_hg38 = function(d, ...) {
  
  geno = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  
  chain = system.file(
    "extdata",
    "hg19ToHg38.over.chain",
    package = "THmisc",
    mustWork = TRUE
  )
  
  lift_data_wrapper(d, chain, geno, ...)
}


