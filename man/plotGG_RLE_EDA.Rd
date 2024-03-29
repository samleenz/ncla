% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plotGG_RLE_EDA}
\alias{plotGG_RLE_EDA}
\title{plotRLE (ggplot edition)}
\usage{
plotGG_RLE_EDA(exprs, group, pal = NULL, title = "")
}
\arguments{
\item{exprs}{row matrix of gene expression (samples as rows)}

\item{group}{vector assigning columns of \code{exprs} to groups to colour by}

\item{pal}{(default = null), character vector of colours the length of
the set of \code{group}. If left default this will be generated by using Set3
which has a max n of 9
@param title a title for the plot}
}
\value{
prints a RLE plot
}
\description{
Make a RLE plot with specific group colourings
}
