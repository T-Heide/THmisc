#' A internal function to print basic file statistics
#'
#' @param f A file.
#'
#' @return NILL
print_file_info = function(f) {
  checkmate::assertFileExists(f, access="r")

  fi = file.info(f)
  cat("\n  File infos:\n")
  cat("   - Name: ", "'", crayon::bold(basename(f)), "'", "\n", sep="")
  cat("   - Created:", crayon::bold(fi$ctime), "\n")
  cat("   - Modified:", crayon::bold(fi$mtime), "\n")
  cat("\n")
}
