#' Internal function to extract VAF information from any 'CollapsedVCF' object.
#'
#' @param d Object of class 'CollapsedVCF'
#'
#' @return A numeric vector containing the VAF data.
get_vaf = function(d) {

  checkmate::assertClass(d, "CollapsedVCF")

  geno_el = names(VariantAnnotation::geno(d)) # elements of the VCF

  if ("VAF" %in% geno_el) {
    # use the existing VAF element
    VAF = VariantAnnotation::geno(d)$VAF

  } else if ("NR" %in% geno_el & "NV" %in% geno_el) {
    # calculate vaf from number of reads and variant reads
    NR = VariantAnnotation::geno(d)$NR; storage.mode(NR) = "numeric"
    NV = VariantAnnotation::geno(d)$NV; storage.mode(NV) = "numeric"
    VAF = NV / NR

  } else {
    #-------------------------------------------------------
    # stop
    #-------------------------------------------------------
    stop("Can't extract VAF data.")
  }

  return(VAF)
}


#' Adds an additional VAF entry to a 'CollapsedVCF' object.
#'
#' @param d A 'CollapsedVCF' object
#'
#' @return A 'CollapsedVCF' object with the VAF added to the genotype information
#' @export
#'
insert_vaf = function(d) {

  checkmate::assertClass(d, "CollapsedVCF")

  if (!"VAF" %in% names(VariantAnnotation::geno(d))) {

    # vaf entry line to add to the header
    vaf_header =
      data.frame(
        Number = ".",
        Type = "Float",
        Description = "Variant allele frequency"
      )

    # generate and add a new header
    old_header = VariantAnnotation::geno(VariantAnnotation::header(d))
    new_header = rbind(old_header, vaf_header)
    rownames(new_header) = c(rownames(old_header), "VAF")
    VariantAnnotation::geno(VariantAnnotation::header(d)) = new_header

    # insert VAF data:
    VariantAnnotation::geno(d)$VAF = get_vaf(d)

  }

  return(d)
}


#' Internal wrapper function to load a Platypus vcf into a format usable for common purposes
#'
#' @param f A Platypus VCF file
#'
#' @return A object of class 'CollapsedVCF'
load_platypus_vcf = function(f) {

  # check arguments
  cl = checkmate::makeAssertCollection()
  checkmate::assertFileExists(f, access="r", add=cl)
  checkmate::reportAssertions(cl)

  # load and modify platpyus vcf file
  d = VariantAnnotation::readVcf(f) %>%
    insert_vaf()


  return(d)
}


#' Function to identify which tool created a given vcf file
#'
#' @param f Path to a existing VCF file.
#'
#' @return A string indentifying the tool that (probably) created the file.
#' @export
indentify_vcf_type = function(f) {

  # check arguments
  checkmate::assertFileExists(f, access="r")

  # load top 10 lines
  d = VariantAnnotation::readVcf(VariantAnnotation::VcfFile(f, yieldSize=10))

  try({ # check if file is a vcf file from platypus
    source = unlist(VariantAnnotation::meta(VariantAnnotation::header(d))$source)
    if (grepl("Platypus_Version_", source))
      return("platypus")
  })

  return("unknown")
}


#' A internal function that summarised the content of a 'CollapsedVCF' object.
#'
#' @param d A 'CollapsedVCF' object.
#'
#' @return NILL
print_vcf_content = function(d) {
  checkmate::assertClass(d, "CollapsedVCF")

  # properties to print
  n_variants = NROW(d)
  n_pass = sum(SummarizedExperiment::rowRanges(d)$FILTER == "PASS")
  frac_pass = round(n_pass/n_variants*100, 1)
  n_chr = length(unique(GenomicRanges::seqnames(SummarizedExperiment::rowRanges(d))))

  # printing
  cat("  VCF object contains:\n")
  cat("   => N=", crayon::bold(n_variants), " variants.\n", sep="")
  cat("   => Of these ", crayon::bold(n_pass, "/", n_variants, sep=""),
      " ", crayon::bold("(", frac_pass,"%)", sep=""), " marked passing.\n", sep="")
  cat("   => On ", crayon::bold(n_chr), " chromosomes\n", sep="")
}


#' Wrapper function to load arbitray VCF files.
#'
#' @param f File name
#' @param verbose Verbosity flag
#'
#' @return A 'CollapsedVCF' object.
#' @export
#'
#' @examples load_vcf_file(system.file("extdata", "platypus.vcf.bgz", package = "THmisc"))
load_vcf_file = function(f, verbose=TRUE) {

  # check arguments
  cl = checkmate::makeAssertCollection()
  checkmate::assertFileExists(f, access="r", add=cl)
  checkmate::assertFlag(verbose, add=cl)
  checkmate::reportAssertions(cl)

  # identify type of vcf file
  if (verbose) {

    # header line
    div = "-"
    div_l = paste0(rep(div, 15), collapse="")
    title = " Loading a VCF file "
    cat(crayon::bold(paste0(div_l, title, div_l, "\n")))
    div_l = paste0(rep(div, 15*2 + nchar(title)), collapse="")

    # print basic file infos
    print_file_info(f)
  }

  if (verbose) cat("  Identifying type... \n")
  vcf_type = indentify_vcf_type(f)
  stopifnot(!is.na(vcf_type))

  # load file
  if (vcf_type == "platypus") {
    if (verbose) cat("   => Type is ", crayon::green("'Platpyus'"), ".\n", sep="")
    data = load_platypus_vcf(f)
  } else if (vcf_type == "mutect2") {
    if (verbose) cat("   => Type is ", crayon::red("Mutect2"), ".\n", sep="")
    cat("\n", crayon::bold(div_l))
    stop("Missing code...\n")
  } else {
    if (verbose) cat("   => Type is ", crayon::red("Unknown"), ".\n", sep="")
    cat("\n", crayon::bold(div_l))
    stop("Unknown VCF file type...\n")  }

  # report content of vcf file:
  if (verbose) {
    cat("\n")
    print_vcf_content(data)
  }

  # print terminal line
  if (verbose)
    cat("\n", crayon::bold(div_l))

  return(data)
}

