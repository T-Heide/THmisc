#' Internal function to extract VAF information from any 'CollapsedVCF' object.
#'
#' @param d Object of class 'CollapsedVCF'
#'
#' @return A numeric vector containing the VAF data.
#' @keywords internal
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

  } else if ("AF" %in% geno_el) {
    # use the existing AF element
    VAF = VariantAnnotation::geno(d)$AF

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



#' Add an additional CNA and AB entries to a 'CollapsedVCF' object.
#'
#' @param mutation_data A 'CollapsedVCF' object
#' @param cna_segments A list of CNA data that can be parsed by the 'get_cnas' function
#'
#' @return A 'CollapsedVCF' object with the CNA and AB added to the genotype information
#' @export
insert_cnas = function(mutation_data, cna_segments) {

  # Add CNA annotation to the header:
  vaf_header =
    data.frame(
      Number=c(".","."),
      Type=c("Integer","String"),
      Description=c("CN","AB")
    )

  old_header = VariantAnnotation::geno(VariantAnnotation::header(mutation_data))
  new_header = rbind(old_header, vaf_header)
  rownames(new_header) = c(rownames(old_header), c("CN","AB"))
  VariantAnnotation::geno(VariantAnnotation::header(mutation_data)) = new_header

  # Insert CNA data:
  VariantAnnotation::geno(mutation_data)$CN =
    get_cnas(mutation_data, cna_segments)

  VariantAnnotation::geno(mutation_data)$AB =
    get_ab_alleles(mutation_data, cna_segments)

  return(mutation_data)
}

#' Adds an additional CCF entry to a 'CollapsedVCF' object.
#'
#' @param mutation_data  A 'CollapsedVCF' object
#' @param purities Ordered purity estimates of all samples in the VCF.
#'
#' @return A 'CollapsedVCF' object with the CCF added to the genotype information
#' @export
insert_ccf = function(mutation_data, purities) {

  stopifnot("VAF" %in% names(VariantAnnotation::geno(mutation_data)))
  stopifnot("CN" %in% names(VariantAnnotation::geno(mutation_data)))
  stopifnot(NCOL(mutation_data) == length(purities))

  # Add CNA annotation to the header:
  added_header = data.frame(Number=".", Type="Float", Description="CCF")
  old_header = VariantAnnotation::geno(VariantAnnotation::header(mutation_data))
  new_header = rbind(old_header, added_header)
  rownames(new_header) = c(rownames(old_header), "CCF")
  VariantAnnotation::geno(VariantAnnotation::header(mutation_data)) = new_header


  # Insert CCF data:
  VariantAnnotation::geno(mutation_data)$CCF =
    VariantAnnotation::geno(mutation_data)$VAF

  for (i in seq_len(NCOL(mutation_data))) {
    vaf = c(unlist(VariantAnnotation::geno(mutation_data)$VAF[,i]))
    cn = c(VariantAnnotation::geno(mutation_data)$CN[,i])
    p = purities[i]
    m = 1 # mutated_alleles

    VariantAnnotation::geno(mutation_data)$CCF[,i] = vaf*(2*(1-p)+cn*p)/(p*m)
  }

  return(mutation_data)
}



#' Internal wrapper function to load a Platypus vcf into a format usable for common purposes
#'
#' @param f A Platypus VCF file
#'
#' @return A object of class 'CollapsedVCF'
#' @keywords internal
load_platypus_vcf = function(f) {

  # check arguments
  cl = checkmate::makeAssertCollection()
  checkmate::assertFileExists(f, access="r", add=cl)
  checkmate::reportAssertions(cl)

  # load and modify platpyus vcf file
  d = VariantAnnotation::readVcf(f) %>%
    split_multiallelic() %>%
    insert_vaf()

  return(d)
}


#' Internal wrapper function to load a Mutect2 vcf into a format usable for common purposes
#'
#' @param f A Platypus VCF file
#'
#' @return A object of class 'CollapsedVCF'
#' @keywords internal
load_mutect2_vcf = function(f) {

  # check arguments
  cl = checkmate::makeAssertCollection()
  checkmate::assertFileExists(f, access="r", add=cl)
  checkmate::reportAssertions(cl)

  # load and modify platpyus vcf file
  d = VariantAnnotation::readVcf(f) %>%
    split_multiallelic() %>%
    insert_vaf()

  return(d)
}


get_csq_parser = function(d) {
  # parse the csq annotation header:
  csq_annot = VariantAnnotation::info(VariantAnnotation::header(d))["CSQ","Description"]
  csq_format = strsplit(gsub("^.*Format: ", "", csq_annot), split="[|]")[[1]]
  (function(l) {return({
    function(e=NULL, f=NULL) {

      if (is.null(e)) {
        cat("Field:\n")
        print(csq_format)
        invisible(csq_format)
      } else {
        if (is.null(f)) {
          stop("Need to pass e and f.\n")
        }
      }

      # if e is a charcter with '|' in each element
      # assume that this is the full csq annoation string
      if (is.character(e)) {
        if (all(grepl("[|]", e))) {
          e = strsplit(e, split="[|]")
        }
      }

      idx = match(f, l)

      if (is.character(e)) {
        res = e[idx]
      } else {
        res = sapply(e, function(e_) e_[idx])
      }

      return(res)
    }
  })})(csq_format)
}


#' Function to split multiallelic sites into multiple variants
#'
#' @param d Object of class 'CollapsedVCF'
#'
#' @return Object of class 'CollapsedVCF' with multiallelic sites split
#' @export
split_multiallelic = function(d) {

  .get_eq_if_default = function(x, y) {
    # modifies alt alleles to match csq output
    # x: ref allele, y: alt allele from header
    if (all(substr(y, 1, 1) == substr(x, 1, 1))) {
      y = substr(y, 2, nchar(y))
      y[y == ""] = "-"
    }

    y
  }

  .split = function(d_cur, i, csq_parser) {

    # info elements:
    el_mult_platypus = c("FR","FS","PP","TR","NF","NR","END")
    el_mult_mutect2 = c("AS_FilterStatus","AS_UNIQ_ALT_READ_COUNT","MPOS","NALOD","POPAF","TLOD")
    for (el in c(el_mult_platypus, el_mult_mutect2)) {
      if (el %in% names(VariantAnnotation::info(d_cur))){
        VariantAnnotation::info(d_cur)[,el] =
          S4Vectors::endoapply(VariantAnnotation::info(d_cur)[,el],  "[", i)
      }
    }

    # info elements per read:
    el_mult_mutect2_r = c("MBQ","MFRL","MMQ","RPA")
    for (el in c(el_mult_mutect2_r)) {
      if (el %in% names(VariantAnnotation::info(d_cur))){
        VariantAnnotation::info(d_cur)[,el] =
          S4Vectors::endoapply(VariantAnnotation::info(d_cur)[,el],"[",c(1,1+i))
      }
    }


    # csq string:
    if ("CSQ" %in% names(VariantAnnotation::info(d_cur))) {

      # get allele infos from vcf
      alts = S4Vectors::lapply(VariantAnnotation::alt(d_cur), as.character)
      ref = as.character(VariantAnnotation::ref(d_cur))
      al_keep = lapply(alts, "[", i)

      # get alt alleles from csq
      csqs = VariantAnnotation::info(d_cur)$CSQ
      al_nums = S4Vectors::lapply(csqs, csq_parser, "ALLELE_NUM")

      if (any(is.na(unlist(al_nums)))) {

        # try modifying the alles so that is matched the ones in CSQ
        # if this fails than CSQ was run with non default options
        al_csq = S4Vectors::lapply(csqs, csq_parser, "Allele")
        exp_csq = mapply(.get_eq_if_default, ref, alts, SIMPLIFY = FALSE)

        if (!all(unlist(mapply("%in%", al_csq, exp_csq)))) {
          # probably the minimal option was used,
          # in this case estimating the csq allele annotation is not
          # done. Complain and exit.
          paste0(      # print warning that ALLELE_NUM are missing and stop.
            "Couldn't find ", sQuote("ALLELE_NUM"), " in the CSQ annotations!\n",
            "  => Please rerun VEP with ", sQuote("--allele_number"), " added."
          ) %>% stop()
        }

        al_nums = mapply(match, al_csq, exp_csq, SIMPLIFY = FALSE)

      }

      # keep the select allele from the csq elements
      VariantAnnotation::info(d_cur)$CSQ =
        mapply("[", csqs, lapply(al_nums, "==", i), SIMPLIFY=0) %>%
        IRanges::CharacterList()

    }


    # gt string:
    if ("GT" %in%  names(VariantAnnotation::geno(d_cur))) {
      VariantAnnotation::geno(d_cur)$GT =
        VariantAnnotation::geno(d_cur)$GT %>%
        gsub(pattern = as.character(i), replacement = "X") %>%
        gsub(pattern = "[1-9]", replacement = "0") %>%
        gsub(pattern = "X", replacement = "1")
    }


    # geno elements:
    for (el in c("NR","NV","AF")) { # only alt
      for (l in seq_len(NCOL(d_cur))){
        if (el %in% names(VariantAnnotation::geno(d_cur))) {
          VariantAnnotation::geno(d_cur)[[el]][,l] =
            S4Vectors::endoapply(VariantAnnotation::geno(d_cur)[[el]][,l],"[",i)
        }
      }
    }

    for (el in c("AD","F1R2","F2R1")) { # ref and alt
      for (l in seq_len(NCOL(d_cur))){
        if (el %in% names(VariantAnnotation::geno(d_cur))) {
          VariantAnnotation::geno(d_cur)[[el]][,l] =
            S4Vectors::endoapply(VariantAnnotation::geno(d_cur)[[el]][,l],"[",c(1,1+i))
        }
      }
    }

    # alt string
    VariantAnnotation::alt(d_cur) =
      S4Vectors::lapply(VariantAnnotation::alt(d_cur),  "[", i) %>%
      Biostrings::DNAStringSetList()

    d_cur
  }

  checkmate::assertClass(d, "CollapsedVCF")

  # keep elements with one allele
  n_elements = S4Vectors::elementNROWS(VariantAnnotation::alt(d))
  wh = n_elements == 1
  d_ret = d[wh,]

  # if no biallelic return data
  n_sites = sum(!wh)
  if (n_sites == 0)
    return(d_ret)

  # elements to drop if not all biallelic
  for (el in c("PL")) {
    if (el %in% names(VariantAnnotation::geno(d))) {
      VariantAnnotation::geno(d_ret)[[el]] = NULL
      VariantAnnotation::geno(d)[[el]] = NULL
    }
  }

  # parse mutli-allelic sites
  d_add = NULL
  csq_parser = get_csq_parser(d)
  max_n_elements = max(n_elements)
  for (i in seq(1, max_n_elements)) {
    idx = n_elements > 1 & n_elements >= i
    if (sum(idx) == 0) next()
    d_cur = .split(d[idx,], i, csq_parser)
    if (is.null(d_add)) {
      d_add = d_cur
    } else {
      d_add = VariantAnnotation::rbind(d_add, d_cur)
    }
  }

  res = VariantAnnotation::rbind(d_ret, d_add)
  res = res[order(res),]

  rownames(res) =
    sprintf(
      "%s:%d_%s/%s",
      GenomicRanges::seqnames(res),
      GenomicRanges::start(res),
      as.character(VariantAnnotation::ref(res)),
      as.character(unlist(VariantAnnotation::alt(res)))
    )

  return(res)
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

  try({ # check if file is a vcf file from mutect2
    meta = VariantAnnotation::meta(VariantAnnotation::header(d))
    version = meta[["MutectVersion"]][,"Value"]
    if (gsub("[.].*", "", version) == "2")
      return("mutect2")
  })

  return("unknown")
}


#' A internal function that summarised the content of a 'CollapsedVCF' object.
#'
#' @param d A 'CollapsedVCF' object.
#'
#' @return NILL
#' @keywords internal
print_vcf_content = function(d) {
  checkmate::assertClass(d, "CollapsedVCF")

  # properties to print
  n_variants = NROW(d)
  n_pass = sum(SummarizedExperiment::rowRanges(d)$FILTER %in% c("PASS","."))
  frac_pass = round(n_pass/n_variants*100, 1)
  n_chr = length(unique(GenomicRanges::seqnames(SummarizedExperiment::rowRanges(d))))

  # printing
  cat("  VCF object contains:\n")
  cat("   => N=", crayon::bold(n_variants), " variants.\n", sep="")
  cat("   => Of these ", crayon::bold(n_pass, "/", n_variants, sep=""),
      " ", crayon::bold("(", frac_pass,"%)", sep=""), " marked passing.\n", sep="")
  cat("   => On ", crayon::bold(n_chr), " chromosomes\n", sep="")
}


#' A internal function to guess the reference genome of a 'CollapsedVCF' object.
#'
#' @param d A 'CollapsedVCF' object.
#' @param verbose Verbosity flag (Default: TRUE)
#'
#' @return Name of the matching reference genome
#' @keywords internal
guess_reference_genome = function(d, verbose=TRUE) {
  checkmate::assertClass(d, "CollapsedVCF")

  # genomes to test
  genomes =
    list(
      "BSgenome.Hsapiens.UCSC.hg38",
      "BSgenome.Hsapiens.UCSC.hg19",
      "BSgenome.Hsapiens.UCSC.hg18"
    )

  if (verbose) {
    cat("  Guessing reference of a vcf file:\n")
  }

  # test if the genomes are available
  avail =
    suppressWarnings(
      sapply(
        genomes,
        require,
        character.only = TRUE,
        quietly = TRUE
      )
    )

  genomes_avail = lapply(genomes[avail], function(x) get(x, envir=asNamespace(x)))

  if (length(genomes_avail) == 0) {
    if (verbose) {
      cat("   => ", crayon::red("No genomes available."), "\n", sep="")
      cat("   => Please install the 'BSgenome.Hsapiens' of the parsed vcf file.\n")
    }
    return(NA)
  }

  if (verbose) {
    n_a = length(genomes_avail)
    n = length(genomes_avail)
    cat("   => Testing ", n_a, "/", n , " available options.\n", sep="")
  }

  # check which genomes match:
  rn = SummarizedExperiment::rowRanges(d)
  rn_ = GenomicRanges::GRanges(
    seqnames = as.character(GenomicRanges::seqnames(rn)),
    ranges = IRanges::IRanges(start(rn), end(rn))
  )

  genomes_matches =
    sapply(genomes_avail, function(g) {
      suppressMessages(suppressWarnings(tryCatch(
        all(as.character(BSgenome::getSeq(g, rn_)) == as.character(rn$REF)),
        error=function(x) return(FALSE)
      )))
    })

  n_m = sum(genomes_matches)
  ref = utils::head(unlist(genomes[avail][genomes_matches]), 1)

  if (verbose) {
    cat("   => Found ", n_m, "/", n_a , " perfectly matching.\n", sep="")
    if (n_m >= 1) {
      cat("   => Reference is ", "'", crayon::green(ref), "'.\n", sep="")
    } else {
      cat(crayon::red("   => Error! Couldn't find reference.\n"), sep="")
    }
  }

  if (length(ref) != 1) {
    return(NA)
  }

  invisible(ref)
}


#' Wrapper function to load arbitrary VCF files.
#'
#' @param f File name
#' @param ... Other optional arguments.
#' @param verbose Verbosity flag (Default: TRUE)
#' @param annot Flag indicating if variant annotations should be added (Default: TRUE)
#'
#' @return A 'CollapsedVCF' object.
#' @export
#'
#' @examples load_vcf_file(system.file("extdata", "platypus.vcf.bgz", package = "THmisc"))
load_vcf_file = function(f, ..., verbose=TRUE, annot=TRUE) {

  # check arguments
  cl = checkmate::makeAssertCollection()
  checkmate::assertFileExists(f, access="r", add=cl)
  checkmate::assertFlag(verbose, add=cl)
  checkmate::assertFlag(annot, add=cl)
  checkmate::reportAssertions(cl)

  # create index if it does not exist
  if (!file.exists(paste0(f, ".tbi")) & grepl("[.]gz$", f)) {
    VariantAnnotation::indexVcf(f)
  }


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

  # indentify the type of the vcf
  if (verbose) cat("  Identifying type... \n")
  vcf_type = indentify_vcf_type(f)
  stopifnot(!is.na(vcf_type))

  if (verbose) {
    if (vcf_type %in% c("unknown")) {
      type = crayon::red(sQuote(Hmisc::capitalize(vcf_type)))
    } else {
      type = crayon::green(sQuote(Hmisc::capitalize(vcf_type)))
    }
    cat("   => Type is ", type, ".\n\n", sep="")
    cat("  Loading data... \n")
  }


  # load data
  if (vcf_type == "platypus") {
    if (verbose) {
      data = load_platypus_vcf(f)
    } else {
      data = suppressMessages(load_platypus_vcf(f))
    }
  } else if (vcf_type == "mutect2") {
    data = load_mutect2_vcf(f)
  } else {
    stop("Unknown VCF file type...\n")
  }


  # drop unused seqleves from data
  seqs_keep = as.character(unique(BSgenome::seqnames(data)))
  GenomeInfoDb::seqlevels(data) = seqs_keep
  GenomeInfoDb::genome(data) = GenomeInfoDb::genome(data)[seqs_keep]


  # report content of vcf file:
  if (verbose) {
    cat("\n")
    print_vcf_content(data)
  }

  # insert annotations
  if (annot) {
    data = add_annot_wrapper(data, verbose=verbose)
  }

  # print terminal line
  if (verbose)
    cat("\n", crayon::bold(div_l), "\n")

  return(data)
}


add_annot_wrapper = function(data, verbose=TRUE) {

  conseq_order = c("frameshift","nonsense","nonsynonymous","synonymous","not translated")

  if ("CSQ" %in% colnames(VariantAnnotation::info(data))) {

    ##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ##  VCF contains VEP annotations    =
    ##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  } else {

    ##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ##  VCF Does not contain annotation =
    ##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # guess reference genome and check that matching TxDb is avail:
    if (verbose) cat("\n")
    rg = guess_reference_genome(data, verbose = verbose)
    if (!is.na(rg)) {
      txdb_n = paste0(gsub("BSgenome[.]", "TxDb.", rg), ".knownGene")
      txdb_a = suppressWarnings(require(txdb_n, character.only = TRUE, quietly = TRUE))
      if (!txdb_a) {
        cat(crayon::red("Missing TxDb package '", "'", sep = ""))
        cat("Please install and try again. Try BiocManager::install('", txdb_n, "')'.\n", sep = "")
        stop("Missing package.")
      }

      # generate annotation object:
      annot = # get annotation
        suppressWarnings(VariantAnnotation::predictCoding(
          query = data,
          subject = get(txdb_n),
          seqSource = get(rg),
        ))

      # print everything missed in the consequence lists
      if (!all(annot$CONSEQUENCE %in% conseq_order)) {
        m_conseq = unique(annot$CONSEQUENCE[!annot$CONSEQUENCE %in% conseq_order])
        stop("Unknown CONSEQUENCES (please update code): ", paste0(m_conseq, collapse = ", "))
      }

      gene_ids_to_name =
        suppressMessages(
          AnnotationDbi::select(
            org.Hs.eg.db::org.Hs.eg.db,
            keys =  unique(annot$GENEID),
            columns = c("ENTREZID", "SYMBOL"),
            keytype = "ENTREZID"
          )
        ) %>% (function(x) magrittr::set_names(x$SYMBOL, x$ENTREZID))


      .get_annot_string = function(mid) {
        annot_c = annot[names(annot) == mid, ,drop=FALSE]

        if (NROW(annot_c) == 0) return(NA)

        if (NROW(annot_c) > 1) {
          if (all(annot_c$CONSEQUENCE == "synonymous")) {
            annot_c = annot_c[1, ]
          } else {
            # keep highest impact
            prior = match(annot_c$CONSEQUENCE, conseq_order)
            annot_c = annot_c[prior == min(prior), ,drop=FALSE]
          }
        }

        if (NROW(annot_c) != 1) {
          warning("Keeping random annotation, need list of cannoical transcripts!\n")
          annot_c = annot_c[1, ,drop=FALSE]
        }

        stopifnot(length(unique(annot_c$REF)) == 1)
        stopifnot(length(unique(unlist(annot_c$ALT))) == 1)

        sprintf(
          "%s p.%s%d%s (%s:%i-%i_%s/%s)",
          gene_ids_to_name[annot_c$GENEID],
          annot_c$REFAA,
          sapply(annot_c$PROTEINLOC, "[[", 1),
          annot_c$VARAA,
          GenomicRanges::seqnames(annot_c),
          GenomicRanges::start(annot_c),
          GenomicRanges::end(annot_c),
          annot_c$REF,
          unlist(annot_c$ALT)
        )

      }

      if (verbose) {
        cat("\n")
        cat("  Generating annotations...\n  ")
      }

      if (verbose & requireNamespace("pbapply", quietly = TRUE)) {
        annot_strings = pbapply::pbsapply(rownames(data), .get_annot_string)
      }  else {
        annot_strings = sapply(rownames(data), .get_annot_string)
      }


      # add annotation of new column to header
      new_el =
        DataFrame(
          Number = 1,
          Type = "Character",
          Description = "Mutation annotation",
          row.names = "ANNOTATION"
        )

      VariantAnnotation::info(VariantAnnotation::header(data)) =
        rbind(VariantAnnotation::info(VariantAnnotation::header(data)), new_el)

      # add mutation annotation to the mutation data
      SummarizedExperiment::rowRanges(data)$ANNOTATION = annot_strings[rownames(data)]

    } else {
      stop("Couldn't identify the BSgenome object to use for annotation.\n")
    }

  }

  return(data)
}



#' Function to drop multiallelic sites from VCFs
#'
#' @param d Object of class 'CollapsedVCF'
#'
#' @return Object of class 'CollapsedVCF' with multiallelic sites split
#' @export
drop_multiallelic = function(d) {
  wh = S4Vectors::elementNROWS(VariantAnnotation::alt(d)) == 1

  if (any(!wh)) {
    warning(sprintf("Dropping %d multiallelic sites!\n", sum(!wh)))
  }

  d[wh,]
}


#' Convert a VCF to a data.frame containing key statistics
#'
#' @param d Object of class 'CollapsedVCF'
#' @param relevant_consequences A character vector defining the relevant consequences to return in the annotations. The possible options are 'high', 'moderate', 'low', 'modifier', 'all' or a combination of these. See https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html for details. (Default: 'high' and 'moderate')
#' @param annotate_all_variants A logical flag indicating if all variants should be included in the annotations or if only the one of the canonical transcript should be keeped. (Default: FALSE)
#' @param annotate_all_genes A logical flag indicating if all genes should be included.  (Default: TRUE)
#' @param include_annot A logical flag indicating if the annotation of the variant should be included. (Default: FALSE)
#' @param trim_canocial A logical flag indicating if the transcript ID of canonical transcripts should be removed. (Default: TRUE)
#' @param genes A character vector of gene names to which the mutation data should be subset. If NULL all genes are included. (Default: NULL)
#'
#' @return A data.frame containing selected information.
#' @export
#'
#' @examples
vcf_to_data_frame = function(d, relevant_consequences=c("high","moderate"), annotate_all_variants=FALSE, annotate_all_genes=TRUE, include_annot=FALSE, trim_canocial=TRUE, genes=NULL) {

  if (is.character(d)) {
    checkmate::assertFileExists(d)
    d = readRDS(d)
  }

  checkmate::assertClass(d, "CollapsedVCF")
  checkmate::assertCharacter(relevant_consequences, null.ok=FALSE, any.missing=FALSE)
  checkmate::assertFlag(annotate_all_variants)
  checkmate::assertFlag(annotate_all_genes)
  checkmate::assertFlag(include_annot)
  checkmate::assertCharacter(genes, null.ok=TRUE, any.missing=FALSE)


  #---------------------------------------------

  # parse selected consequences
  rel_csq_in = relevant_consequences
  relevant_consequences = c()

  for (i in seq_along(rel_csq_in)) {
    if (rel_csq_in[i] %in% names(THmisc:::csq_consequence_mapper)) {
      added_csqs = THmisc:::csq_consequence_mapper[[rel_csq_in[i]]]
    } else {
      if (rel_csq_in[i] %in% unlist(THmisc:::csq_consequence_mapper)) {
        added_csqs = rel_csq_in[i]
      } else {
        alt_csq = gsub(" ", "_", tolower(rel_csq_in[i]))
        if (alt_csq %in% unlist(THmisc:::csq_consequence_mapper)) {
          added_csqs = alt_csq
        } else{
          msg = paste0(
            "Unkown consequence ",
            sQuote(rel_csq_in[i]),
            ".\n",
            "Use ",
            sQuote("THmisc:::csq_consequence_mapper"),
            " to see valid ones\n."
          )
          stop(msg)
        }
      }
    }
    relevant_consequences = c(relevant_consequences, added_csqs)
  }

  #---------------------------------------------

  empty_result =
    data.frame(
      mutation=character(),
      sample=character(),
      variant=character(),
      CN=numeric(),
      CCF=numeric(),
      VAF=numeric(),
      NV=numeric(),
      NR=numeric(),
      AB=character(),
      type=factor(levels = c("SNV","MNV","InDel")),
      gene=character()
    )

  if (NROW(d) == 0)
    return(empty_result)

  #---------------------------------------------

  has_annot = "ANNOTATION" %in% colnames(S4Vectors::elementMetadata(SummarizedExperiment::rowRanges(d)))
  has_csq = "CSQ" %in% names(VariantAnnotation::info(d))

  #---------------------------------------------

  if (has_csq) {

    .get_feature = THmisc:::get_csq_parser(d)

    .replace_aa_codes = function(v) {
      p = c("Ter", seqinr::aaa())
      r = c("*", seqinr::a())
      for (i in seq_along(p))
        v = gsub(p[i], r[i], v)
      return(v)
    }

    .gen_annot = function(x, trim_canocial) {

      if (length(x) == 0)
        return("")

      paste(sapply(x, function(b) {
        gene = .get_feature(b, "SYMBOL")

        conseq = .replace_aa_codes(.get_feature(b, "HGVSp"))

        if (trim_canocial & conseq != "")
          conseq = strsplit(conseq , ":")[[1]][2]

        if (conseq == "") {# no protein code, try coding
          conseq = .get_feature(b, "HGVSc")

          if (trim_canocial & conseq != "")
            conseq = strsplit(conseq , ":")[[1]][2]

          if (conseq != "")
            conseq = paste0(conseq, ", ")

          conseq = paste0(conseq, stringr::str_to_title(gsub("_", " ", .get_feature(b, "Consequence"))))
        }

        conseq = gsub("%3", "", conseq)

        paste0(gene, " (", conseq, ")")
      }), collapse="; ")
    }

    .get_most_relevant_csq = function(x) {

      if (annotate_all_variants) {

        wh_annotate = rep(TRUE, length(x))

      } else {

        conseq = strsplit(.get_feature(x, "Consequence"), "&")
        is_trans = .get_feature(x, "Feature_type") == "Transcript"
        is_canon = .get_feature(x, "Feature") %in% THmisc:::canonical_transcripts$transcript

        max_ = function(x) { if (all(is.na(x))) return(NA); max(x, na.rm=TRUE) }
        most_relevant = sapply(lapply(conseq, match, relevant_consequences), max_)
        most_relevant[is.infinite(most_relevant)] = NA
        if (all(is.na(most_relevant))) {
          is_relevant = rep(FALSE, length(most_relevant))
        } else {
          is_relevant = most_relevant == max(most_relevant, na.rm=TRUE)
          is_relevant = sapply(is_relevant, isTRUE)
        }

        # identify elements to use for annotation
        wh_annotate = is_relevant & is_canon & is_trans
        trim_canocial = TRUE
        if (sum(wh_annotate) == 0) {
          wh_annotate = is_relevant & is_trans
          trim_canocial = FALSE
          if (sum(wh_annotate) == 0) {
            wh_annotate = is_relevant
          }
          if (!annotate_all_genes & any(wh_annotate)){
            wh_annotate = which(wh_annotate)[1]
          }
        }

      }

      annot = .gen_annot(x[wh_annotate], trim_canocial)

      #}, error=function(e){print(e); NA})

    }

  }

  #---------------------------------------------

  if (!is.null(genes)) {

    if (has_csq) {

      .keep_el =  function(x) any(unlist(x) %in% genes)

      # drop irrelevant lines
      csq_sym = lapply(VariantAnnotation::info(d)$CSQ, .get_feature, "SYMBOL")
      csq_sym_spl = csq_sym %>% lapply(strsplit, "[%]")
      lines_keep = sapply(csq_sym_spl, .keep_el)
      d = d[lines_keep,]

      # drop irrelevant annotations
      csq_sym_spl = csq_sym_spl[lines_keep]
      csq_keep = lapply(csq_sym_spl, sapply, .keep_el)
      VariantAnnotation::info(d)$CSQ =
        mapply(function(x, y) x[y], VariantAnnotation::info(d)$CSQ, csq_keep, SIMPLIFY = FALSE) %>%
        CharacterList()

    } else if (has_annot) {

      gene = gsub(" [(].*", "", SummarizedExperiment::rowRanges(d)$ANNOTATION)
      wh_keep = sapply(gene, function(x) any(x %in% genes))
      d = d[wh_keep,]
      warning("Does not contain CSQ. Will not filter consequences.\n")

    } else {
      stop("Can't filter for genes since CSQ strings are missing. Add code?\n")
    }
  }

  if (NROW(d) == 0)
    return(empty_result)

  #---------------------------------------------

  if (has_csq & include_annot) {
    mut_annot = sapply(VariantAnnotation::info(d)$CSQ, .get_most_relevant_csq)
    names(mut_annot) = rownames(d)
  } else if (has_annot) {
    mut_annot = SummarizedExperiment::rowRanges(d)$ANNOTATION
    names(mut_annot) = rownames(d)
  } else {
    mut_annot = rep(NA, NROW(d))
    names(mut_annot) = rownames(d)
  }

  .retabulate = function(x) {
    x_ = unlist(x)
    dim(x_) = dim(x)
    x_ = data.frame(x_)
    dimnames(x_) = dimnames(x)
    cbind(data.frame(id=rownames(x_)), x_)
  }

  to_merge = list()
  for (e in c("CN","CCF","VAF","NV","NR","AB")) {
    if (e %in% names(VariantAnnotation::geno(d))) {
      to_merge[[e]] =
        .retabulate(VariantAnnotation::geno(d)[[e]]) %>%
        reshape::melt(id.var="id")
    } else {
      if (e == "NR" & "DP" %in% names(VariantAnnotation::geno(d))) {
        to_merge[[e]] =
          .retabulate(VariantAnnotation::geno(d)$DP) %>%
          reshape::melt(id.var="id")
      } else if (e == "NV" & "AD" %in% names(VariantAnnotation::geno(d))) {
        nv = apply(VariantAnnotation::geno(d)$AD, 2, sapply, "[", 2)
        dim(nv) = dim(d)
        dimnames(nv) = dimnames(d)
        to_merge[[e]] = .retabulate(nv) %>%  reshape::melt()
      } else if (e == "VAF" & "AF" %in% names(VariantAnnotation::geno(d))) {
        af = VariantAnnotation::geno(d)$AF
        to_merge[[e]] = .retabulate(af) %>%  reshape::melt(id.vars="id")
      } else if (e %in% c("AB","CN","CCF")) { # cna data are optional
        na_res = matrix(NA, nrow = NROW(d), ncol = NCOL(d))
        dimnames(na_res) = dimnames(d)
        to_merge[[e]] = reshape::melt(na_res)
      }
    }
  }

  for (i in seq_along(to_merge)[-1]) {
    checkmate::assertTRUE(all(to_merge[[1]][,1:2] == to_merge[[i]][,1:2]))
    to_merge[[i]] = to_merge[[i]][,-c(1:2)]
  }

  # merge CN, VAF and AB etc. counts data:
  .get_gene = function(x) {
    as.character(x) %>%
      strsplit(split="[; ]",fixed = TRUE) %>%
      lapply(gsub, pattern=" [(][&.+->?_,: A-Za-z0-9*]+[)]", replacement="") %>%
      lapply(gsub, pattern=" [pc][.].*", replacement="") %>%
      lapply(gsub, pattern=" ", replacement="") %>%
      sapply(paste0, collapse=";")
  }

  data_bound =
    do.call(what=cbind, to_merge) %>%
    magrittr::set_colnames(c("mutation","sample", names(to_merge))) %>%
    dplyr::mutate(variant = mut_annot[as.character(mutation)]) %>%
    dplyr::mutate(gene = .get_gene(variant)) %>%  # set mutation types:
    dplyr::mutate(type = THmisc::get_mutation_type(mutation))  # set mutation types:


  checkmate::assertSetEqual(colnames(empty_result), colnames(data_bound))

  #---------------------------------------------

  # order AB counts (keep up to cn of 20):
  n_max = max(c(4, max(nchar(as.character(data_bound$AB)), na.rm=TRUE)))

  ab_levels =
    unlist(lapply(0:n_max, function(x)
      lapply(seq(from = 0, to = floor(x / 2)), function(y)
        paste0(c(rep("A", x - y), rep("B", y)), collapse = ""))))

  ab_levels = ab_levels[order(ab_levels)]
  ab_levels = ab_levels[order(nchar(ab_levels))]

  stopifnot(all(as.character(data_bound$AB) %in% ab_levels | is.na(data_bound$AB)))
  data_bound$AB = factor(data_bound$AB, ab_levels, ordered=TRUE)


  return(data_bound[,colnames(empty_result)])
}


get_annotation_from_csq = function(d) {

  # parse the variant annotations:
  get_feature = get_csq_parser(d)
  annot_csq = lapply(VariantAnnotation::info(d)$CSQ, strsplit, split="[|]")


  # find relevant annotations:
  annot_csq_relevant = sapply(annot_csq, function(a) {

    conseq = strsplit(sapply(a, get_feature, "Consequence"), "&")

    is_trans = sapply(a, get_feature, "Feature_type") == "Transcript"
    is_canon = sapply(a, get_feature, "Feature") %in% canonical_list$transcript
    is_relevant = sapply(lapply(conseq, "%in%", relevant_consequences), any)

    wh_annotate = is_trans & is_canon & is_relevant
    trim_canocial = TRUE
    if (sum(wh_annotate) == 0) {
      wh_annotate = is_trans & is_relevant
      trim_canocial = FALSE
      if (sum(wh_annotate) == 0) {
        wh_annotate = is_relevant
      }
      if (!annotate_all & any(wh_annotate)){
        wh_annotate = which(wh_annotate)[1]
      }
    }

    annots_genes =
      paste(sapply(a[wh_annotate], function(b) {
        gene = get_feature(b, "SYMBOL")
        conseq = replace_aa_codes(get_feature(b, "HGVSp"))
        if (conseq == "")
          conseq = get_feature(b, "HGVSc")
        if (trim_canocial)
          conseq = strsplit(conseq , ":")[[1]][2]
        conseq = gsub("%3", "", conseq)
        paste0(gene, " (", conseq, ")")
      }), collapse="; ")

    if (length(annots_genes) == 0) {
      return("")
    } else {
      return(annots_genes)
    }
  })


}

replace_aa_codes = function(v) {
  p = c("Ter", seqinr::aaa())
  r = c("*", seqinr::a())
  for (i in seq_along(p))
    v = gsub(p[i], r[i], v)
  return(v)
}

#' Plot a histogram of mutations contained in a VCF across CN states.
#'
#' @param d Object of class 'CollapsedVCF'
#' @param what String containing 'VAF' or 'CCF' that defines what to plot. (Default: 'VAF')
#' @param purity A named numeric vector of sample purities. Used to calculate expected cluster positions and to create labels. (Default: NULL)
#' @param ploidy  A named numeric vector of sample ploidies. Used to calculate expected cluster positions and to create labels. (Default: NULL)
#'
#' @return A ggplot object containing mutation histograms
#' @export
#' @import ggplot2
#'
plot_vaf_distribution = function(d, what="VAF", purity=NULL, ploidy=NULL) {

  checkmate::assertDataFrame(d)
  checkmate::assertTRUE(all(c("mutation","sample","CN","VAF","CCF","AB","type") %in% colnames(d)))
  checkmate::assertCharacter(what, max.len = 1, min.len = 1)
  checkmate::assertChoice(what, c("VAF","CCF"))

  checkmate::assertNumeric(purity, null.ok = TRUE, names = "unique")
  checkmate::assertNumeric(ploidy, null.ok = TRUE, names = "unique")

  # plot data for each sample
  d$value = d[,what]
  plot_list = list()

  for (csample in unique(as.character(d$sample))) {

    title = csample
    subtitle = ""

    if (csample %in% names(purity)) {
      subtitle = paste0(subtitle, "Purity: ", round(purity[csample] * 100), "%")
    }

    if (csample %in% names(ploidy)) {
      if (subtitle != "") subtitle = paste0(subtitle, ", ")
      subtitle = paste0(subtitle, "Ploidy: ", round(ploidy[csample]))
    }


    cn_labels = labeller(
      CN = function(s) {
        s = if_else(as.numeric(s) > 5, "CN > 4", paste("CN =", as.character(s)))
      }
    )

    ab_labels = labeller(
      AB = function(string) {
        paste0("CN: ", as.character(string), "")
      }
    )


    data_plot=
      data.frame(d) %>%
      dplyr::filter(sample == csample) %>%
      dplyr::filter(value > 0) %>%
      dplyr::mutate(AB = factor(AB)) %>%
      mutate(AB=factor(AB, levels(AB)[nchar(levels(AB)) <= 4], ordered=TRUE)) %>%
      dplyr::filter(!is.na(AB))


    if (NROW(data_plot) == 0) {
      next()
    }

    mutations_per_cn =
      table(data_plot$AB) %>%
      reshape2::melt() %>%
      magrittr::set_colnames(c("AB","N")) %>%
      mutate(N=paste0("N = ", N))


    if (csample %in% names(purity)) {

      ab_states =
        mutations_per_cn %>%
        dplyr::filter(N != "N = 0") %>%
        dplyr::select(AB) %>% unlist()


      .get_peak_pos = function(ab) {

        p = purity[csample]
        ab = as.character(ab)
        cn = nchar(ab)
        n_states = sum(strsplit(ab, "")[[1]] == "A")
        am = seq_len(n_states)

        if (what == "VAF") {
          data.frame(
            AB = rep(ab, n_states),
            value = am * p / (2 * (1 - p) + cn * p),
            af = am,
            row.names = NULL
          )
        } else if (what == "CCF") {
          data.frame(
            AB = rep(ab, n_states),
            value = am,
            af = am,
            row.names = NULL
          )
        } else {
          NULL
        }
      }

      peak_position_per_cn = do.call(rbind, lapply(ab_states, .get_peak_pos))
    } else {
      peak_position_per_cn = NULL
    }

    plot =
      data_plot %>%
      ggplot(aes(x=value)) +
      geom_histogram(aes(fill=type), bins=50,  alpha=1.0)

    if (!is.null(peak_position_per_cn)) {
      plot = plot +
        geom_vline(data=peak_position_per_cn, aes(xintercept=value, group=af), linetype=2, color="gray30")
    }

    plot = plot +
      geom_text(data=mutations_per_cn, aes(label=N, x=1, y=Inf), hjust=1, vjust=2, color="black") +
      facet_wrap(~AB, scales="free", drop=FALSE, labeller=ab_labels) +
      xlim(0, 1) +
      theme(plot.subtitle=element_text(hjust=0.5)) +
      ggtitle(title, subtitle) +
      scale_fill_brewer(palette="Set1", drop=FALSE) +
      theme(panel.spacing = unit(1, "lines")) +
      xlab(what) +
      ylab("Number") +
      labs(fill="Mutation type")


    plot_list[[csample]] = ggplotGrob(plot)

  }

  return(plot_list)
}


#' Plot CN states of somatic sides across samples.
#'
#' @param d Object of class 'CollapsedVCF'
#' @param min_freq Minimum relative size to include a group of CN states. (Default: 0.01)
#'
#' @return A ggplot object containing mutation histograms
#' @export
#' @import ggplot2
#'
plot_cn_states_across_samples = function(data, min_freq = 0.01) {

  cn_stat_mat = with(data, tapply(as.character(AB), list(mutation, sample), c))

  # drop samples with no CN infos
  wh_drop = apply(is.na(cn_stat_mat), 2, all)
  cn_stat_mat = cn_stat_mat[,!wh_drop]

  # convert to string
  cn_stats_v = apply(cn_stat_mat, 1, paste0, collapse=" ")
  cn_stats_tab = table(cn_stats_v)

  # summary data frame
  cn_states =
    data.frame(
      cn_state = names(cn_stats_tab),
      freq = as.numeric(cn_stats_tab / sum(cn_stats_tab)),
      count = as.numeric(cn_stats_tab)
    )

  plot =
    cn_states %>%
    dplyr::filter(freq >= min_freq) %>%
    ggplot(aes(y=cn_state, x="N", fill=freq, label=count)) +
    geom_tile() +
    geom_text() +
    scale_fill_distiller(palette=8, direction=1) +
    xlab("") +
    ylab("CN state")

  return(plot)
}


#' Plot a heatmap showing VAF of specific mutations across samples.
#'
#' @param d Data frame containing the data to plot.
#' Has to include the columns 'sample', 'NV', 'NR', 'VAF', 'gene', 'variant'.
#' Can optionally include the column 'label' which is used to label tiles (defaults to 'NV'/'NR').
#' Can optionally include the column 'group' which is used to split samples into tiles.
#' @param max_label_width A integer defining the maxmimum length of labels. Labels longer than this are trimmed.
#'
#' @rdname plot_driver_vaf_heatmap
#' @return
#' @export
#' @import ggplot2
#'
plot_driver_vaf_heatmap = function(d, max_label_width = 50) {

  checkmate::assertDataFrame(d)
  checkmate::assertTRUE(all(c("sample","NV","NR","VAF","gene","variant") %in% colnames(d)))

  if ("group" %in% colnames(d)) {
    .all_equal = function(x) length(unique(x)) <= 1
    unique_groups_per_sample = all(sapply(split(d$group, d$sample), .all_equal))
    checkmate::assertTRUE(unique_groups_per_sample)
  }

  if (NROW(d) == 0)
    return(NULL)

  #---------------------------------------------

  gene_spl = lapply(strsplit(d$gene, "[&,; ]"), function(x) unique(x[x!=""]))
  d$gene = sapply(gene_spl, paste0, collapse=", ")

  #---------------------------------------------

  # Trim excessively long labels
  d$variant = as.character(d$variant)
  char_length = nchar(d$variant)
  wh_to_long = char_length > max_label_width
  if (any(wh_to_long)) {
    trimmed_label = strtrim(d$variant[wh_to_long], width = max_label_width - 4)
    d$variant[wh_to_long] = paste(trimmed_label, "...")
  }

  #---------------------------------------------

  # Make fixed width axis labels (pad with spaces)
  char_length = nchar(d$variant)
  if (!any(char_length == max_label_width)) {
    wh_max = char_length == max(char_length)
    append = paste0(rep(" ", max_label_width - max(char_length)), collapse = "")
    d$variant[wh_max] = paste0(append, d$variant[wh_max])
  }

  #---------------------------------------------

  if (!"label" %in% colnames(d)) {
    d$label = paste0(d$NV, "/", d$NR)
  }

  #---------------------------------------------

  plot_of_drivers =
    ggplot(d, aes(x = sample, y = variant, fill = VAF, label=label)) +
    cowplot::theme_cowplot() +
    geom_bin2d(size = 0.2, color = "gray20") +
    geom_text(size = 2.5, alpha = 0.8) +
    facet_grid(gene ~ ., scales = "free", space = "free") +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
    background_grid() +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    xlab("Sample") +
    ylab("Mutation")

  # optionally split panel by column 'group'
  if ("group" %in% colnames(d)) {
    plot_of_drivers = plot_of_drivers +
      facet_grid(gene ~ group, scales = "free", space = "free")
  }

  return(plot_of_drivers)
}



#' Wrapper around plot_driver_vaf_heatmap function to save the image with resonable dimensions
#'
#' @param d Data frame containing the data to plot.
#' @param out_file Name of the image file to create.
#' @param genes A character vector of genes to subset data to. (Default: NULL)
#' @param ... Arguments passed to plot_driver_vaf_heatmap. See ?plot_driver_vaf_heatmap.
#'
#' @rdname plot_driver_vaf_heatmap
#' @return
#' @export
#'
save_driver_vaf_heatmap = function(d, out_file, genes=NULL, ...) {

  checkmate::assertDataFrame(d)
  checkmate::assertTRUE(all(c("sample","NV","NR","VAF","gene","variant") %in% colnames(d)))
  checkmate::assertCharacter(genes, null.ok = TRUE)

  if (!is.null(genes)) {
    spl = strsplit(d$gene, "[&,; ]")
    keep = sapply(spl, function(x) isTRUE(any(unlist(x) %in% genes)))
    d = d[keep,]
  }

  plot = plot_driver_vaf_heatmap(d, ...)

  # find good width
  n_sample_sets = length(unique(plot$data$group))
  n_samples = length(unique(plot$data$sample))
  lab_width = max(nchar(plot$data$variant))
  width =  0.13 * (n_sample_sets - 1) + 0.15 * n_samples + 0.127 * lab_width + 1.95

  # find good height
  n_genes = length(unique(plot$data$gene))
  n_mutations = length(unique(paste0(plot$data$variant, "...", plot$data$gene)))
  height = 0.13 * (n_genes - 1) + 0.17 * n_mutations + 0.27 + 2.43

  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  ggsave(out_file, plot, width=width, height=height, limitsize=FALSE)

  invisible(plot)
}


#' Create a plot showing if a gene is mutated in a given sample
#'
#' @param d Data frame containing the data to plot.
#' @param crit_value Critical value of column 'use' to consider a sample to be mutated (i.e., d[,use] > crit_value)
#' @param genes # Genes to include in the plot.
#' @param use Column of the data frame containing information if a gene is mutated.
#' @param order_alphabetic Flag indicating if rows and columns should be ordered alphabetically.
#' @param gene_to_gr # Named character vector or labeller function used to convert gene ids to a group label.
#' @param keep_all_genes # Flag indicating if genes unmutaed in all samples should be included in the plot.
#' @param xlab # Label for the x-axis (e.g., the genes)
#' @param order_patients # Predefined order of samples, if logical TRUE then the calculated order is returned instead.
#'
#' @return
#' @export
#'
plot_gene_mut_heatmap = function(d, crit_value=0.1, genes=NULL, use="CCF", order_alphabetic=FALSE, gene_to_gr=NULL, keep_all_genes=FALSE, xlab="", order_patients=NULL) {

  checkmate::assertDataFrame(d)
  checkmate::assertTRUE(all(c("patient","sample",use,"gene","variant") %in% colnames(d)))
  checkmate::assertCharacter(genes, null.ok = TRUE)
  checkmate::assertNumeric(d[,use], null.ok = FALSE, any.missing = FALSE)
  checkmate::assertFlag(order_alphabetic)
  if (!isTRUE(order_patients))
    checkmate::assertCharacter(order_patients, null.ok=TRUE)

  #---------------------------------------------

  if (is.character(gene_to_gr)) {
    gene_to_gr = (function(m) {force(m); return(function(x) m[x])})(gene_to_gr)
  } else {
    if (!is.function(gene_to_gr) & !is.null(gene_to_gr)) {
      stop("'gene_to_gr' should be a map or labeller function.")
    }
  }

  #---------------------------------------------

  if (!"patient_gr" %in% colnames(d)) {
    d$patient_gr = ""
  }

  #---------------------------------------------

  gene_spl = lapply(strsplit(d$gene, "[&,; ]"), function(x) unique(x[x!=""]))
  d$gene = sapply(gene_spl, paste0, collapse=", ")
  d$value = d[,use]
  d$patient_gr = factor(d$patient_gr)

  #---------------------------------------------

  if (!is.null(genes)) { # subset to specific genes
    spl = strsplit(d$gene, "[&,; ]")
    keep = sapply(spl, function(x) isTRUE(any(unlist(x) %in% genes)))
    d = d[keep,]

    if (keep_all_genes) {
      d$gene = factor(d$gene, unique(genes, unique(d$gene)), ordered = TRUE)
    }
  }

  #---------------------------------------------

  if (NROW(d) == 0)
    return(NULL)

  #---------------------------------------------

  # pad patient ids with spaces depending on tissue
  tissues = unique(d$patient_gr)
  mt = match(d$patient_gr, tissues)
  pad = sapply(mt, function(n) paste0(rep(" ", n-1), collapse = ""))
  d$patient = paste0(pad, d$patient)

  #---------------------------------------------

  # combinations to tabulate
  genes = if (is.factor(d$gene)) levels(d$gene) else unique(d$gene)
  patients = unique(d$patient)
  tissues = unique(d$patient_gr)
  samples_per_set = lapply(split(as.character(d$sample), list(d$patient, d$patient_gr)), unique)
  tissues_per_pat = lapply(split(as.character(d$patient_gr), d$patient), unique)

  # data frame for tabulation
  tab_data = expand.grid(gene=genes, patient=patients, tissue=tissues, stringsAsFactors = FALSE)
  for (i in seq_along(tissues_per_pat)) { # drop non existing tissue patient combinations
    pat = names(tissues_per_pat)[i]
    val_t = tissues_per_pat[[i]]
    wh_drop = tab_data$patient == pat & !tab_data$tissue %in% val_t
    tab_data = tab_data[!wh_drop,]
  }
  tab_data$fraction_mutated = NA
  tab_data$mean_value = NA
  tab_data$status = NA
  tab_data$type = NA

  for (i in seq_len(NROW(tab_data))) {

    # get view of input data
    wh_c = d$gene == as.character(tab_data$gene[i]) &
      d$patient == as.character(tab_data$patient[i]) &
      d$patient_gr == as.character(tab_data$tissue[i])

    d_c = d[wh_c,]
    if (NROW(d_c) == 0) next()

    # mapper from variant to mutation
    wh_ndup = !duplicated(d_c$mutation)
    mut_to_var = d_c$variant[wh_ndup] %>%
      make.unique() %>%
      magrittr::set_names(d_c$mutation[wh_ndup])

    # expand variant data to matrix
    grouping = list(as.character(d_c$sample), as.character(d_c$mutation))
    value_matrix = tapply(d_c$value, grouping, unique)

    # assert that all samples are present
    exp_smp = with(tab_data[i, ], samples_per_set[[paste0(patient, ".", tissue)]])
    stopifnot(all(exp_smp %in% rownames(value_matrix)))

    # calculate values to tabulate
    frac_mut = apply(value_matrix > crit_value, 2, mean)
    mean_value = apply(value_matrix, 2, mean)

    # other propeties in order of above table
    mt = match(gsub("[.][0-9]+$", "", colnames(value_matrix)), d_c$mutation)
    types = d_c$type[mt]

    # selection of most relevant mutation
    wh_use = frac_mut == max(frac_mut)
    if (sum(wh_use) > 1 & is.ordered(types)) {
        wh_use = wh_use & max(types[wh_use], na.rm=TRUE) == types
    }

    if (sum(wh_use) > 1) {
      wh_use = wh_use & mean_value == max(mean_value[wh_use])
    }

    if (sum(wh_use) > 1) { # take care not to call sample if only one is selected
      wh_use = sample(which(wh_use), size = 1)
    }

    # save data
    tab_data$fraction_mutated[i] = frac_mut[wh_use]
    tab_data$mean_value[i] = mean_value[wh_use]
    tab_data$status[i] = if_else(frac_mut[wh_use] == 1, "clonal", "sub-clonal")
    tab_data$type[i] = as.character(types[wh_use])
  }

  if (!keep_all_genes) {
    wh = tab_data$gene %in% tab_data$gene[which(tab_data$fraction_mutated > 0)]
    tab_data = tab_data[wh,]
  }

  #---------------------------------------------
  .is_mut = function(x) if (all(is.na(x))) return(0) else return(max(x))
  mut_mat = with(tab_data, tapply(fraction_mutated, list(patient, gene), .is_mut))

  if (order_alphabetic) {
    order_genes = sort(as.character(tab_data$gene), decreasing = TRUE)
    order_pat = sort(as.character(tab_data$patient), decreasing = FALSE)
  } else {
    order_genes = names(sort(apply(mut_mat>0, 2, sum), decreasing = TRUE))
    ord_list = split(t(mut_mat[,order_genes]), seq_len(NCOL(mut_mat)))
    order_pat = rownames(mut_mat)[do.call(order, c(ord_list, decreasing=TRUE))]
  }

  if (!is.null(order_patients)) {

    if (isTRUE(order_patients))
      return(order_pat)

    order_pat = order_patients
  }

  tab_data$patient = factor(tab_data$patient, order_pat, ordered = TRUE)
  tab_data$gene = factor(tab_data$gene, rev(order_genes), ordered = TRUE)

  #---------------------------------------------

  if (!is.null(gene_to_gr)) {
    tab_data$gene_gr = gene_to_gr(as.character(tab_data$gene))
  }

  #---------------------------------------------

  n_gene_mutated = apply(mut_mat > 0 , 2, sum)
  n_patients = NCOL(mut_mat)
  label_mut = sprintf("n = %d", n_gene_mutated)
  y_pos_mut = match(names(n_gene_mutated), order_genes)
  mutation_labs = data.frame(y=y_pos_mut, x=n_patients+1, label=label_mut)

  #---------------------------------------------

  tab_data$plot_value = tab_data$fraction_mutated

  #---------------------------------------------

  lf = labeller(.default = Hmisc::capitalize)

  plot_heatmap =
    tab_data %>%
    dplyr::mutate(plot_value = ifelse(plot_value > 0, plot_value, NA)) %>%
    dplyr::filter(is.na(plot_value)) %>%
    ggplot(aes(x=patient, y=gene, fill=plot_value)) +
    theme_cowplot() +
    geom_tile(color="gray0", size=0.2) +
    scale_fill_continuous(na.value="white", guide="none") +
    facet_grid(.~as.character(tissue), scales="free", space="free", drop=FALSE, labeller = lf) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
    theme(legend.position="bottom", legend.box="vertical") +
    theme(strip.text = element_text(size=11), strip.background = element_blank()) +
    xlab(xlab) +
    ylab("") +
    labs(color="")

  #---------------------------------------------

  if (is.factor(d$type)) {
    uniq_types = levels(d$type)
  } else {
    uniq_types = na.omit(sort(unique(d$type)))
  }
  fill_colors = c("#bd0026","#08589e", "#810f7c", rainbow(n=10))
  stopifnot(length(uniq_types) <= length(fill_colors))

  fill_guide = guide_legend(
    title.position = "top",
    title.hjust = 0.5,
    override.aes = list(shape = NA)
  )

  for (i in rev(seq_along(uniq_types))) {
    plot_heatmap =
      plot_heatmap +
      ggnewscale::new_scale_fill() +
      geom_tile(
        data = dplyr::filter(tab_data, type == uniq_types[i] & plot_value > 0),
        aes(fill = plot_value),
        color = "gray0",
        size = 0.2
      ) +
      scale_fill_gradient(
        na.value = "white",
        low = "white",
        high = fill_colors[i],
        limits = c(-0.2, 1.1),
        guide = fill_guide
      ) +
      labs(fill=paste0("Fraction of samples mutated (", uniq_types[i], ")"))
  }

  #---------------------------------------------

  plot_heatmap =
    plot_heatmap +
    geom_point(
      data = dplyr::filter(tab_data, status == "clonal"),
      aes(color = "Clonal"),
      size = 0.8
    ) +
    scale_color_manual(
      values = c("Clonal" = "darkorange")
    ) +
    guides(
      color =
        guide_legend(
          title.position = "top",
          title.hjust = 0.5
        )
    ) +
    labs(color=" ")

  #---------------------------------------------

  plot_heatmap =
    plot_heatmap + # ensure correct order of labels
    scale_x_discrete(breaks = order_pat,  labels = gsub(" ","", order_pat))
  #  scale_y_discrete(limits = rev(order_genes))

  #---------------------------------------------

  if ("gene_gr" %in% colnames(tab_data)) {
    # if there is a gene group in the data split plot by it

    plot_heatmap =
      plot_heatmap +
      facet_grid(
        gene_gr ~ tissue,
        scales = "free",
        space = "free",
        drop = FALSE,
        labeller = lf
      )
  }

  #---------------------------------------------

  gt = ggplotGrob(plot_heatmap)
  for (i in which(grepl("strip-[tr]", gt$layout$name))){
    gt$grobs[[i]]$layout$clip = "off"
  }

  return(gt)
}


#' Adds a 'cluster' id to a mutation data frame
#'
#' @param d A mutation data frame.
#' @param ccf_cutoff Cutoff of the CCF to consider sample mutated [default:0.25].
#' @param small_frac Relative size of the smallest cluster to label [default:0.01].
#'
#' @return A mutation data frame with a cluster id column 'clust_id' added.
#' @export
add_cluster_infos = function(d, ccf_cutoff=0.25, small_frac=0.01) {

  checkmate::assert_data_frame(d)
  checkmate::assertSubset(c("CCF", "mutation","sample"), colnames(d))
  checkmate::assertNumber(ccf_cutoff, lower = 0, upper = Inf)
  checkmate::assertNumber(small_frac, lower = 0, upper = 1)

  # spread mut to matrix
  mut_mat = tapply(d$CCF > ccf_cutoff, list(d$mutation, d$sample), as.numeric)
  mut_mat = mut_mat[,!apply(is.na(mut_mat), 2, all)]

  # order case ids
  mut_mat = mut_mat[,order(apply(mut_mat, 2, sum, na.rm=TRUE))]
  smp_ids = unique(c(colnames(mut_mat), as.character(d$sample)))
  d$sample = factor(d$sample, smp_ids)

  # collapse to id
  mut_id = apply(mut_mat, 1, paste0, collapse="")
  n_large = max(c(0, length(mut_id) * small_frac))
  large_groups = names(which(table(mut_id) > n_large))
  mut_id[!mut_id %in% large_groups] = "small"

  # add as factor to data
  ord_mut_id = unique(mut_id[order(apply(mut_mat, 1, sum), decreasing=FALSE)])
  ord_mut_id = rev(unique(c("small", rev(ord_mut_id))))
  d$mut_id = factor(mut_id[d$mutation], ord_mut_id, ordered = TRUE)

  cl_order = paste0("C", seq_along(levels(d$mut_id))-1)
  cl_order[rev(levels(d$mut_id)) == "small"] = "S"

  d$clust_id = factor(d$mut_id, levels(d$mut_id), rev(cl_order), ordered=TRUE)

  return(d)
}


#' Function to plot mutation data as a heatmap
#'
#' @param d A mutation data frame.
#' @param value Name of numeric id to plot. (default: "VAF")
#' @param annot Optional annotation of mutations. Names have to be mutation ids and the values the annotation.
#' @param title Optional title to add to plot.
#'
#' @return A ggplot object.
#' @export
#' @import ggplot2
#'
plot_mutation_heatmap = function(d, value="VAF", annot=NULL, title="") {

  checkmate::assert_data_frame(d)
  checkmate::assertSubset(c(value, "mutation", "sample", "clust_id"), colnames(d))
  checkmate::assertCharacter(annot, names = "named", null.ok = TRUE)
  checkmate::assertString(title)

  d$value_plot = d[,value]

  pl =
    dplyr::filter(d, !is.na(sample)) %>%
    ggplot(aes(y=mutation, x=sample, fill=value_plot)) +
    facet_grid(clust_id~., scale="free", space="free") +
    geom_tile() +
    xlab("Sample") +
    ylab("") +
    labs(fill=value) +
    scale_fill_gradient(
      low = "white",
      high = "darkred",
      limits = c(0, 1),
      oob = scales::squish,
      na.value = "gray99"
    ) +
    theme(
      axis.text.x = element_text(angle=90, vjust=0.5),
      axis.text.y = element_blank(),
      strip.text.y = element_text(size=10, angle=0),
      axis.line.y=element_blank()
    ) +
    ggtitle(
      title, subtitle=sprintf("N = %d", length(unique(d$mutation)))
    )


  # split by group variable
  if ("group" %in% colnames(d))
    pl = pl + facet_grid(clust_id~group, scale="free", space="free")


  # add mutation annotations
  if (!is.null(annot)) {
    annot = annot[names(annot) %in% d$mutation]

    pl = pl + scale_y_discrete(breaks=names(annot), labels=annot) +
      theme(axis.text.y = element_text(size=10, color=alpha("black", 0.8)))
  }

  return(pl)
}


#' Function to plot mutation data as a histogram
#'
#' @param d A mutation data frame.
#' @param value Name of numeric id to plot [default: "VAF"].
#' @param max_value Maximum value fo the value column. [default: Inf].
#' @param title Optional title to add to plot.
#'
#' @return A ggplot object.
#' @export
#' @import ggplot2
#'
plot_mutation_histogram = function(d, value="VAF", max_value=Inf, title="") {

  checkmate::assert_data_frame(d)
  checkmate::assertSubset(c(value, "mutation", "sample"), colnames(d))
  checkmate::assertString(title)
  checkmate::assertNumber(max_value, na.ok = FALSE, null.ok = FALSE)

  d$value_plot = d[,value]
  if (!"clust_id" %in% colnames(d)) d$clust_id = NA

  n_groups = length(na.omit(unique(d$clust_id)))
  col_vals = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(n_groups)

  # Histograms
  pl = d %>%
    mutate(value_plot=if_else(value_plot>max_value, max_value, value_plot)) %>%
    ggplot(aes(x=value_plot, fill=clust_id, group=paste0(sample, clust_id))) +
      geom_histogram(binwidth=0.025, alpha=1, position="identity") +
      facet_wrap(~as.character(sample), ncol=4, scales = "free_y") +
      scale_fill_manual(values=col_vals, na.value = "gray30") +
      xlim(0.01, NA) +
      xlab(value) +
      ylab("Count") +
      labs(fill="Cluster")+
      ggtitle(ggtitle(title, sprintf("N = %d", length(unique(d$mutation)))))

  return(pl)
}



#' Function to plot intermutation distances
#'
#' @param d A mutation data frame. Can be obtained using 'vcf_to_data_frame' from a 'CollapsedVCF' object.
#' @param title Optional title to add to plot. . (Default: "")
#' @param geno A 'BSgenome' object that can be used to set limits and order of chromosomes. (Default: NULL)
#'
#' @return A ggplot object.
#' @export
#' @import ggplot2
#'
plot_inter_mutation_distance = function(d, title="", geno=NULL) {

  checkmate::assertTRUE(inherits(d, c("data.frame", "CollapsedVCF")))
  checkmate::assertString(title)

  if (inherits(d, "data.frame")) {
    checkmate::assertTRUE("mutation" %in% colnames(d))
    mut_ids = unique(d$mutation)
  } else if (inherits(d, "CollapsedVCF")) {
    mut_ids = d
  } else {
    stop()
  }

  mut_data = THmisc::get_mutation_df(mut_ids)
  mut_data$type = THmisc::get_mutation_type(mut_ids)
  mut_data$transition = THmisc:::get_transition_type(mut_ids)
  n_muts = NROW(mut_ids)

  # insert inter mutation distance
  get_mm_dist = function(x) sapply(seq_along(x), function(i) min(abs(x[i]-x[c(i-1,i+1)]), na.rm=TRUE))
  starts_per_chr = split(mut_data$start, mut_data$chr)
  mut_data$mut_mut_dist = unlist(lapply(starts_per_chr, get_mm_dist))

  if (!is.null(geno)) {

    # order of chromosomes
    wh = GenomicRanges::seqnames(geno) %in% mut_data$chr
    chr_ord = GenomicRanges::seqnames(geno)[wh]
    mut_data$chr = factor(mut_data$chr, chr_ord, ordered=TRUE)

    # mark of end of each chromosome
    end_marks =
      data.frame(
        chr = chr_ord,
        start = seqlengths(geno)[chr_ord],
        end = NA,
        ref = NA,
        alt = NA,
        type = "SNV",
        mut_mut_dist = max(mut_data$mut_mut_dist),
        transition = NA
      )

    mut_data =
      rbind(
        mut_data,
        end_marks
      )

  }

  # Plot waterfall plot
  pl =
    mut_data %>%
    dplyr::filter(type == "SNV") %>%
    ggplot(aes(x=start, y=mut_mut_dist, color=transition)) +
    geom_point(alpha=0.8, size=0.5) +
    facet_grid(~chr, scales="free_x", space="free_x") +
    scale_y_log10() +
    xlab("Position") +
    labs(color="Type") +
    theme(axis.text.x=element_blank(), axis.ticks=element_blank()) +
    theme(strip.text.x=element_text(angle=90)) +
    background_grid() +
    scale_color_manual(values = trans_colors, breaks=names(trans_colors)) +
    ylab("Inter-mutation distance") +
    guides(colour=guide_legend(override.aes=list(size=1, alpha=1))) +
    ggtitle(title, sprintf("N = %d", n_muts))

  return(pl)
}
