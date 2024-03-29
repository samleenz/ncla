% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plotGG_UMAP_eda}
\alias{plotGG_UMAP_eda}
\title{Plot UMAP (ggplot edition)}
\usage{
plotGG_UMAP_eda(dat, group, pal = NULL, title = "")
}
\arguments{
\item{dat}{umap object, output of \code{umap::umap()}}

\item{group}{vector assigning samples from \code{dat} to groups to colour by}

\item{pal}{(default = null), character vector of colours the length of
the set of \code{group}. If left default this will be generated by using Set1
which has a max n of 9}

\item{title}{title of the plot}
}
\description{
Plot UMAP (ggplot edition)
}
