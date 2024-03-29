% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cn_analysis.R
\name{calc_avg_cn_per_bin}
\alias{calc_avg_cn_per_bin}
\title{Calculate average copy-number values in bins.}
\usage{
calc_avg_cn_per_bin(cn_data, windows)
}
\arguments{
\item{cn_data}{A list of data frames or a single data frame that contains interval definitions (start/start.pos, end/end.pos, chr/chromosome columns) and copy-number values (CN/CNt).}

\item{windows}{A data.frame that can be converted to a GRanges object defining the window intervals.}
}
\value{
A data.frame containing all columns of the 'windows' argument and average copy-number values of all entries in 'cn_data'.
}
\description{
Calculate average copy-number values in bins.
}
