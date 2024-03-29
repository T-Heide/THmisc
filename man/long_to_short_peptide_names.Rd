% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snvs.R
\name{long_to_short_peptide_names}
\alias{long_to_short_peptide_names}
\title{Convert long peptide names to short ones.}
\usage{
long_to_short_peptide_names(annot)
}
\arguments{
\item{annot}{A character vector of annotations with peptide abbreviations}
}
\value{
Also a character vector, put with short peptide names
}
\description{
Convert long peptide names to short ones.
}
\examples{
long_to_short_peptide_names(c("TerAlaVal"))
}
