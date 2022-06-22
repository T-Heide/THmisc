#' Get the mutation type of mutations
#'
#' @param d A CollapsedVCF object or a mutation string.
#'
#' @return A named factor vector containing the mutation type (SNV, MNV or InDel).
#' @export
#'
#' @examples get_mutation_type("chr1:757771_A/G")          # A single-nucleotide variant (SNV).
#' @examples get_mutation_type("1:757771_A/G")             # Also a SNV.
#' @examples get_mutation_type("chr1:757771_AT/AG")        # A multinucleotide variant.
#' @examples get_mutation_type("chr1:757771-757772_AT/AG") # Also a MNV.
#' @examples get_mutation_type("chr1:757771_ATG/A")        # A InDel.
#' @examples get_mutation_type(c("chr1:757770_A/G","chr1:757771_ATG/A")) # Two mutations
get_mutation_type = function(d) {

  mut_df = get_mutation_df(d)

  with(mut_df, {
    # mutation types
    type = rep("SNV", NROW(d))               # default is SNV
    type[nchar(ref) != 1] = "MNV"            # MNV are more than 1 nucleotide long,
    type[nchar(ref) != nchar(alt)] = "InDel" # but not an InDel

    # set names and order levels:
    names(type) = rownames(mut_df)

    type = factor(type, c("SNV","MNV","InDel"), ordered=TRUE)
    return(type)
  })
}


#' Get the transition type of mutations
#'
#' @param d A CollapsedVCF object or a mutation string.
#'
#' @return A named factor vector containing the mutation type (SNV, MNV or InDel).
#' @export
#'
#' @examples get_mutation_type("chr1:757771_A/G")          # A single-nucleotide variant (SNV).
#' @examples get_mutation_type("1:757771_A/G")             # Also a SNV.
#' @examples get_mutation_type("chr1:757771_AT/AG")        # A multinucleotide variant.
#' @examples get_mutation_type("chr1:757771-757772_AT/AG") # Also a MNV.
#' @examples get_mutation_type("chr1:757771_ATG/A")        # A InDel.
#' @examples get_mutation_type(c("chr1:757770_A/G","chr1:757771_ATG/A")) # Two mutations
get_transition_type = function(d) {

  mut_df = get_mutation_df(d)
  mut_df$type = get_mutation_type(mut_df)

  with(mut_df, {

    rev_comp = ref %in% c("G", "A")
    ref[rev_comp] = c("A"="T", "C"="G", "G"="C", "T"="A")[ref[rev_comp]]
    alt[rev_comp] = c("A"="T", "C"="G", "G"="C", "T"="A")[alt[rev_comp]]

    trans_type = paste0(ref, ">", alt)
    trans_type[type != "SNV"] = NA

    trans_types = c("C>A","C>G","C>T","T>A","T>C","T>G")
    return(factor(trans_type, trans_types, ordered=TRUE))
  })
}


#' Get a data frame of mutation from various input formats
#'
#' @param d A CollapsedVCF object or a mutation string.
#'
#' @return A data.frame containing mutation informations (chromosome, start, end and ref. & alt. alleles).
#' @export
#'
#' @examples get_mutation_df("chr1:757771_A/G")          # A single-nucleotide variant (SNV).
#' @examples get_mutation_df("1:757771_A/G")             # Also a SNV.
#' @examples get_mutation_df("chr1:757771_AT/AG")        # A multinucleotide variant.
#' @examples get_mutation_df("chr1:757771-757772_AT/AG") # Also a MNV.
#' @examples get_mutation_df("chr1:757771_ATG/A")        # A InDel.
#' @examples get_mutation_df(c("chr1:757770_A/G","chr1:757771_ATG/A")) # Two mutations
get_mutation_df = function(d) {

  if (is.data.frame(d)) {
    if (all(c("chr","start","end","ref","alt") %in% colnames(d))) {
      return(d)
    }
  }

  if (is.factor(d)) {
    d = as.character(d)
  }

  if (is.character(d)) {

    # check input is a valid mutation string
    mut_regex = "^(chr)?[0-9a-zA-z]+:[0-9]+(-[0-9]+)?_[ACGTN]+/[ACGTN]+$"
    checkmate::assertCharacter(d)
    checkmate::assertTRUE(all(grepl(mut_regex, d)))

    # get ref and alt alleles
    d_spl = strsplit(d, "[:_/]")
    chr = sapply(d_spl, "[", 1)
    pos = sapply(d_spl, "[", 2)
    start = as.numeric(gsub(".*-", "", pos))
    end = as.numeric(gsub("-.*", "", pos))
    ref = sapply(d_spl, "[", 3)
    alt = sapply(d_spl, "[", 4)
    mut_ids = d

  } else if (inherits(d, "CollapsedVCF")){

    # check input only contains one alt per row
    enr  = S4Vectors::elementNROWS(VariantAnnotation::alt(d))
    checkmate::assertTRUE(all(enr == 1))

    # get ref and alt alleles
    chr = GenomicRanges::seqnames(d)
    start = GenomicRanges::start(d)
    end = GenomicRanges::end(d)
    ref = VariantAnnotation::ref(d)
    alt = unlist(VariantAnnotation::alt(d))
    mut_ids = rownames(d)

  } else {
    stop("Please pass a 'CollapsedVCF' object or a mutation string as input.\n")
  }

  df =
    data.frame(
      chr = chr,
      start = start,
      end = end,
      ref = ref,
      alt = alt,
      row.names = NULL
    )

  if (!any(duplicated(mut_ids))) {
    rownames(df) = mut_ids
  }

  return(df)

}


#' Complemental strand of a DNA sequence with preserved direction.
#'
#' @param x DNA base sequence.
#' @return Complemental strand with identical orientation (3'->5' and 3'->5').
#' @examples baseRevComplement("ACCCG")
#' @keywords internal
#' @export
baseRevComplement =
  function(x){

    # Return NULL if x is NULL:
    if (is.null(x)) return(NULL)

    # Convert x to vector of characters:
    x = as.vector(as.character(x))

    # Test if all elements of x containe a valid sequences:
    if (!all(grepl("^([ATCG]*)$", x)))
      stop("Elements of x must only contain DNA bases [ACGT].")

    x = chartr(x, old = "ATCG", new = "TAGC")
    out = strsplit(x, split = "") %>% lapply(rev) %>% sapply(paste, collapse = "")

    names(out) = NULL
    return(out)
  }


createMutTypes = function(c_str, a_str) {

  mutation_types =
    c(
      "A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C",
      "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T",
      "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "A[C>G]A", "A[C>G]C",
      "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T",
      "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C",
      "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
      "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C",
      "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
      "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C",
      "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T",
      "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "A[T>C]A", "A[T>C]C",
      "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T",
      "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C",
      "T[T>C]G", "T[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
      "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C",
      "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"
      )

  c = Biostrings::DNAStringSet(c_str)
  a = Biostrings::DNAStringSet(a_str)

  width = nchar(c_str)
  p_center = (width-1)/2+1

  wh = substr(c, p_center, p_center) %in% c("G", "A")
  c[wh] = Biostrings::reverseComplement(c[wh])
  a[wh] = Biostrings::reverseComplement(a[wh])

  mt = paste0(
    substr(c, 1, p_center - 1),
    "[",
    substr(c, p_center, p_center),
    ">",
    as.character(a),
    "]",
    substr(c, p_center + 1, width)
  )

  if (all(width == 3)) {
    mt = factor(mt, mutation_types, ordered=TRUE)
  }

  return(mt)
}


#' Create detailed mutation annotation.
#'
#' @param x A CollapsedVCF object or a mutation string.
#' @param geno A BSgenome object containing the reference genome.
#' @param n_bases_context Number of base-pairs to use as context.

#' @return A data.frame containing detailed mutation annotation (chromosome, start, end and ref. & alt. alleles, mutation type, transition + context that can be used for signature analysis and transition type).
#' @export
#'
#' @examples if (require("BSgenome.Hsapiens.UCSC.hg38")) {
#'   THmisc:::annotate_mutation_types(
#'     c("chr1:11925339_A/G", "chr1:11925592_C/T"),
#'     BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'   )
#' }
annotate_mutation_types = function(x, geno, n_bases_context=1) {

  # basic annotation
  muts = get_mutation_df(x)
  muts$mutation_type = get_mutation_type(muts)

  #
  wh_snv = muts$mutation_type == "SNV"
  muts_gr = as(muts[wh_snv,c("chr","start","end")], "GRanges")
  cont = as.character(BSgenome::getSeq(geno, muts_gr + n_bases_context))

  # check reported reference allele
  width = nchar(cont)
  p_center = (width-1)/2+1
  ref = substr(cont, p_center, p_center)
  if (any(ref != muts$ref[wh_snv]))
    stop("Mismatch in reported reference and geno. Wrong genome defined via 'geno'?")

  # add substitiution type
  sub_types = createMutTypes(cont, muts$alt[wh_snv])
  muts$mutation_context = factor(NA, levels=levels(sub_types))
  muts$mutation_context[wh_snv] = sub_types
  muts$substitution = gsub("[]].*", "", gsub(".*[[]", "", muts$mutation_context))

  return(muts)
}


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

