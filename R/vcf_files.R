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



get_csq_parser = function(d) {
  # parse the csq annotation header:
  csq_annot = info(header(d))["CSQ","Description"]
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
      
      if (is.character(e)) {
        e = strsplit(e, split="[|]")
      }
      
      idx = match(f, l)
      res = sapply(e, function(e_) e_[idx])
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
  
  .split = function(d_cur, i) {
    
    # info elements:
    for (el in c("FR","PP","TR","NF","NR","FS")) {
      if (el %in% names(VariantAnnotation::info(d_cur))){
        VariantAnnotation::info(d_cur)[,el] = 
          S4Vectors::endoapply(VariantAnnotation::info(d_cur)[,el],  "[", i)
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
        CharacterList()
      
    }
    
    
    # gt string:
    VariantAnnotation::geno(d_cur)$GT = 
      VariantAnnotation::geno(d_cur)$GT %>% 
      gsub(pattern = as.character(i), replacement = "X") %>% 
      gsub(pattern = "[1-9]", replacement = "0") %>% 
      gsub(pattern = "X", replacement = "1")
    
    
    for (l in seq_len(NCOL(d_cur))){
      VariantAnnotation::geno(d_cur)$NR[,l] = 
        S4Vectors::endoapply(VariantAnnotation::geno(d_cur)$NR[,l],  "[", i)
      
      VariantAnnotation::geno(d_cur)$NV[,l] = 
        S4Vectors::endoapply(VariantAnnotation::geno(d_cur)$NV[,l],  "[", i)
    }
    
    VariantAnnotation::alt(d_cur) = # alt reads
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

  # parse mutli-allelic sites  
  d_add = NULL
  csq_parser = get_csq_parser(d)
  max_n_elements = max(n_elements)
  for (i in seq(1, max_n_elements)) {
    idx = n_elements > 1 & n_elements >= i
    if (sum(idx) == 0) next()
    d_cur = .split(d[idx,], i)
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

  genomes_avail = lapply(genomes[avail], function(x) get(x))

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

  genomes_matches =
    sapply(genomes_avail, function(g) {
      suppressMessages(suppressWarnings(tryCatch(
        all(as.character(BSgenome::getSeq(g, rn)) == rn$REF),
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
  if (!file.exists(paste0(f, ".tbi"))) {
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

  if (verbose) cat("  Identifying type... \n")
  vcf_type = indentify_vcf_type(f)
  stopifnot(!is.na(vcf_type))

  # load file
  if (vcf_type == "platypus") {
    if (verbose) {
      cat("   => Type is ", crayon::green("'Platpyus'"), ".\n\n", sep="")
      cat("  Loading data... \n")
      data = load_platypus_vcf(f)
    } else {
      data = suppressMessages(load_platypus_vcf(f))
    }
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

  conseq_order = c("frameshift","nonsense","nonsynonymous","synonymous")

  if (FALSE) {

    ##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ##  VCF contains VEP annotations    =
    ##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    stop()

  } else {

    ##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ##  VCF Does not contain annotation =
    ##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # guess reference genome and check that matching TxDb is avail:
    if (verbose) cat("\n")
    rg = guess_reference_genome(data, verbose = verbose)
    if (!is.na(rg)) {
      txdb_n = paste0(gsub("BSgenome[.]", "TxDb.", rg), ".knownGene")
      txdb_a = suppressWarnings(require(txdb_n, character.only = TRUE, quietly =
                                          TRUE))
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
        stop("Unknown CONSEQUENCES: ", paste0(m_conseq, collapse = ", "))
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
        annot_c = annot[names(annot) == mid, ]

        if (NROW(annot_c) == 0) return(NA)

        if (NROW(annot_c) > 1) {
          if (all(annot_c$CONSEQUENCE == "synonymous")) {
            annot_c = annot_c[1, ]
          } else {
            # keep highest impact
            prior = match(annot_c$CONSEQUENCE, conseq_order)
            annot_c = annot_c[prior == min(prior), ]
          }
        }

        if (NROW(annot_c) != 1) {
          warning("Keeping random annotation, need list of cannoical transcripts!\n")
          annot_c = annot_c[1, ]
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
      SummarizedExperiment::rowRanges(data)$ANNOTATION = annot_strings[rownames(data)]

    } else {
      stop("Couldn't identify the BSgenome object to use for annotation.\n")
    }

  }

  return(data)
}
