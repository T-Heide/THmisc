% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snvs.R
\name{annotate_mutation_types}
\alias{annotate_mutation_types}
\title{Create detailed mutation annotation.}
\usage{
annotate_mutation_types(x, geno, n_bases_context = 1)
}
\arguments{
\item{x}{A CollapsedVCF object or a mutation string.}

\item{geno}{A BSgenome object containing the reference genome.}

\item{n_bases_context}{Number of base-pairs to use as context.}
}
\value{
A data.frame containing detailed mutation annotation (chromosome, start, end and ref. & alt. alleles, mutation type, transition + context that can be used for signature analysis and transition type).
}
\description{
Create detailed mutation annotation.
}
\examples{
if (require("BSgenome.Hsapiens.UCSC.hg38")) {
  THmisc:::annotate_mutation_types(
    c("chr1:11925339_A/G", "chr1:11925592_C/T"),
    BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  )
}
}
