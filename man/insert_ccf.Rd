% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_files.R
\name{insert_ccf}
\alias{insert_ccf}
\title{Adds an additional CCF entry to a 'CollapsedVCF' object.}
\usage{
insert_ccf(mutation_data, purities)
}
\arguments{
\item{mutation_data}{A 'CollapsedVCF' object}

\item{purities}{Ordered purity estimates of all samples in the VCF.}
}
\value{
A 'CollapsedVCF' object with the CCF added to the genotype information
}
\description{
Adds an additional CCF entry to a 'CollapsedVCF' object.
}
