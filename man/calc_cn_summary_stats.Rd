% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cn_analysis.R
\name{calc_cn_summary_stats}
\alias{calc_cn_summary_stats}
\title{Calculate copy-number alteration summary stats}
\usage{
calc_cn_summary_stats(d, window_infos = NULL, ...)
}
\arguments{
\item{d}{A data.frame containing CN data. Must contain a 'sample' or 'sample_barcode', a 'patient' , 'CN' or 'value' (the copy-number value), and 'name' (the CN window name).}

\item{window_infos}{A optional interval annotation data frame of windows. Must contain "chr","start","end","name","window" columns.}

\item{...}{arguements passed to `findAmplificationBins`}
}
\value{
A list of  data.frame containing summary statics about how often a genomic window ('name') is lost, gained or amplified in each 'patient'.
}
\description{
Calculate copy-number alteration summary stats
}
