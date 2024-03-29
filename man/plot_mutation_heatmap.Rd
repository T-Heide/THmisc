% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_files.R
\name{plot_mutation_heatmap}
\alias{plot_mutation_heatmap}
\title{Function to plot mutation data as a heatmap}
\usage{
plot_mutation_heatmap(d, value = "VAF", annot = NULL, title = "")
}
\arguments{
\item{d}{A mutation data frame.}

\item{value}{Name of numeric id to plot. (default: "VAF")}

\item{annot}{Optional annotation of mutations. Names have to be mutation ids and the values the annotation.}

\item{title}{Optional title to add to plot.}
}
\value{
A ggplot object.
}
\description{
Function to plot mutation data as a heatmap
}
