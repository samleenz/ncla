% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing.R
\name{transformData}
\alias{transformData}
\title{transformData: apply a transformtation to row matrix}
\usage{
transformData(x, f)
}
\arguments{
\item{x}{the row matrix}

\item{f}{transformation to apply, one of:
\itemize{
\item none: return the same matrix
\item log2: log-2 transform with pseudo-count of one
\item scale: scale and center
\item pcaN  : PCA where N PCs are kept
}}
}
\value{
the transformed data matrix
}
\description{
transformData: apply a transformtation to row matrix
}
