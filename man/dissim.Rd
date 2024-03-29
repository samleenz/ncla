% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dissimilarity.R
\name{dissim}
\alias{dissim}
\title{dissim}
\usage{
dissim(x, f, g = NULL)
}
\arguments{
\item{x}{row matrix}

\item{f}{dissimilarity measure to use; one of
\itemize{
\item "L1"
\item "L2"
\item "cosine"
\item "pearson"
\item "spearman"
\item "JSD"
\item "bhattacharyya"
\item "TSD"
note: jensen-shannon and bhattacharyya distance use relative
frequencies to calculate probability vectors and will
fail if data contains negative values
}}

\item{g}{optional, the gene space to use for TSD, if NULL -
all genes will be used.}
}
\value{
n x n dissimilarity matrix of x using function type f
}
\description{
Calculate a dissimilarity matrix using a specified metric (or measure)
}
