% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_files.R
\name{plot_mutation_histogram}
\alias{plot_mutation_histogram}
\title{Function to plot mutation data as a histogram}
\usage{
plot_mutation_histogram(d, value = "VAF", max_value = Inf, title = "")
}
\arguments{
\item{d}{A mutation data frame.}

\item{value}{Name of numeric id to plot [default: "VAF"].}

\item{max_value}{Maximum value fo the value column. [default: Inf].}

\item{title}{Optional title to add to plot.}
}
\value{
A ggplot object.
}
\description{
Function to plot mutation data as a histogram
}
