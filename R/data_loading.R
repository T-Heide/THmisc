#' Convience function to load Sequenza fit files.
#'
#' @param files Vector of files to load.
#' @param verbose Verbosity flag (Default: FALSE).
#'
#' @return Returns a data.frame containing per sample summary stats (purity, ploidy and average ploidy).
#' @export
load_sequenza_purity_files = function(files, verbose=FALSE) {
  
  # check arguments
  checkmate::assertFileExists(files, access = "r")
  checkmate::assertFlag(verbose)
  
  
  # arguments for loader function
  args = list(
    X = files,
    FUN = readr::read_tsv,
    col_types = "ddd",
    skip = 2,
    n_max = 1,
    col_names = c("purity","ploidy","ploidy_mean"), 
    progress=FALSE
  )
  
  # load data files
  if (verbose & requireNamespace("pbapply", quietly = TRUE)) {
    data = do.call(what=pbapply::pblapply, args) 
  } else {
    data = do.call(what=lapply, args)
  }
  data = do.call(what=rbind, data)
  
  if (!is.null(names(files))) {
    data = cbind(sample=names(files), data)
  } else {
    data = cbind(sample=gsub("_confints_CP[.]txt$", "", basename(files)), data)
  }
  
  
  # return results
  invisible(data)
}


#' Convience function to load Sequenza segment files.
#'
#' @param files Vector of files to load.
#' @param verbose Verbosity flag (Default: FALSE).
#'
#' @return Returns a list containing per sample data.
#' @export
load_sequenza_segment_files = function(files, verbose = FALSE) {
  
  # check arguments
  checkmate::assertFileExists(files, access = "r")
  checkmate::assertFlag(verbose)
  
  
  # load data files
  col_types = "ciiddddidiiid"
  if (verbose & requireNamespace("pbapply", quietly = TRUE)) {
    data = pbapply::pblapply(files, readr::read_tsv, col_types=col_types, progress=FALSE)
  } else {
    data = lapply(files, readr::read_tsv, col_types=col_types, progress=FALSE)
  }
  
  if (is.null(names(files))) {
    names(data) = gsub("_segments[.]txt$", "", basename(files))
  }
  
  # return results
  invisible(data)
}
