#' Calculate RCD
#'
#' Given a matrix x with samples in rows and genes in columns calculate the
#' Ranking Correlation Distance defined my Manatakis et al. 2021
#'
#' Optionally give a subset of genes to use as the "atlas genes" as vector g
#'
#' @param x
#' @param g
#'
#' @return
#' @export
#'
#' @examples
calcRCD <- function(x, g=NULL){

  # calculate sample (row-wise) ranks
  #   breaks ties with first
  #   highest expression should have max rank so use -x
  x_rank <- apply(x, 1, function(v) rank(-v, ties.method = "first"))

  # subset for the genes in G
  if(is.null(g)){
    x_rank_g <- x_rank
  } else {
    x_rank_g <- x_rank[g, ]
  }

  # calculate pairwise RCCs (pearson correlation of Ranks)
  RCC <- cor(x_rank_g, method = "pearson")

  # calcualte RCD
  #   sqrt(1-ReLU(RCC))
  RCD <- sqrt(1 - pmax(RCC, 0))

  return(RCD)
}

#' Calculate srJSD
#'
#' Calculate square-root of the Jensen-Shannon divergence
#' for rows of matrix x
#'
#' Optionally give a subset of genes to use as the "atlas genes" as vector g
#'
#' @param x
#' @param g
#'
#' @return
#' @export
#'
#' @examples
calcSRJSD <- function(x, g = NULL){
  # subset genes of x by the atlas genes if given
  if(! is.null(g)){
    x <- x[, g]
  }

  jsd <- suppressMessages(philentropy::JSD(x, est.prob = "empirical"))
  srjsd <- sqrt(jsd)

  # catch for if only two observations are being compared,
  # make the 2x2 distance matrix instead of a single distance
  if(length(srjsd) == 1){
    srjsd <- matrix(c(0, srjsd, srjsd, 0), nrow = 2)
  }
  dimnames(srjsd) <- list(rownames(x), rownames(x))

  return(srjsd)
}

#' Transcriptome Signature distance
#'
#' Given a gene expression matrix, x, with samples in rows calculate the
#' pairwise Transcriptome Signature Distances as per Manatakis et al. 2021
#'
#' Optionally give a subset of genes to use as the "atlas genes" as vector g
#'
#' @param x
#' @param g
#'
#' @return
#' @export
#'
#' @examples
calcTSD <- function(x, g = NULL){
  TSD <- (0.5 * calcSRJSD(x,g)) + (0.5 * calcRCD(x,g))

  return(TSD)
}
