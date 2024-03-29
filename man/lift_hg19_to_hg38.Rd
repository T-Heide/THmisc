% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lifting_between_references.R
\name{lift_hg19_to_hg38}
\alias{lift_hg19_to_hg38}
\title{Function to lift mutations or objects between genome versions}
\usage{
lift_hg19_to_hg38(d, ...)
}
\arguments{
\item{d}{Input data to work on. Options are i) a list of mutation ids, ii) a GRanges object.}

\item{...}{Additional arguments passed to lift_data_wrapper (see ?lift_data_wrapper for details)}
}
\value{
A object with a structure similar to d with updated position data.
}
\description{
Function to lift mutations or objects between genome versions
}
\examples{
lift_hg19_to_hg38("chr2:148654004_T/C")
}
