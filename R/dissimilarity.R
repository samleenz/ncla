#' dissim
#'
#' Calculate a dissimilarity matrix using a specified metric (or measure)
#'
#' @param x row matrix
#' @param f dissimilarity measure to use; one of
#'   - "L1"
#'   - "L2"
#'   - "cosine"
#'   - "pearson"
#'   - "spearman"
#'   - "JSD"
#'   - "bhattacharyya"
#'   - "TSD"
#'   note: jensen-shannon and bhattacharyya distance use relative
#'   frequencies to calculate probability vectors and will
#'   fail if data contains negative values
#' @param g optional, the gene space to use for TSD, if NULL -
#'   all genes will be used.
#'
#' @return n x n dissimilarity matrix of x using function type f
#' @export
#'
#' @examples
dissim <- function(x, f, g = NULL){
  if(f == "L1"){
    d <- philentropy::distance(x, "manhattan", mute.message = TRUE)
  } else if(f == "L2"){
    d <- philentropy::distance(x, "euclidean", mute.message = TRUE)
  } else if(f == "cosine"){
    d <- acos(cosine(x)) / pi
  } else if(f == "pearson"){
    d <- (1 - stats::cor(t(x), method = 'pearson')) / 2
  } else if(f == "spearman"){
    d <- (1 - stats::cor(t(x), method = 'spearman')) / 2
  } else if(f == "jensen-shannon"){
    d <- philentropy::distance(
      x,
      "jensen-shannon",
      est.prob = "empirical",
      mute.message = TRUE
      )
  } else if(f == "bhattacharyya"){
    d <- philentropy::distance(
      x,
      "bhattacharyya",
      est.prob = "empirical",
      mute.message = TRUE
      )
  } else if(f == "F5"){ # fractional minkowksi distnace, p = 0.5
    d <- philentropy::distance(x, "minkowski", p = 0.5, mute.message = TRUE)
  } else if(f == "TSD"){# Transcriptome Signature Distance
    d <- calcTSD(x, g)
  } else {
    stop("f is not a valid option")
  }

  # for philentropy functions
  # if matrix passed to philentropy::distance() is 2 row a single value
  # will be returned instead of a 2x2 distance matrix
  if(length(d) == 1){
    matrix(c(0, d, d, 0), nrow = 2)
  } else{
    d
  }
}


#' cosine: calculate the cosine similarity of rows of matrix x
#' function form modeled on coss() from
#' https://stackoverflow.com/questions/45377700/fast-computation-of-106-cosine-vector-similarities-in-r
#' diag(z) <- 1 added as this is a similarity matrix by definitation
#' if diag not set to 1 cannot calc acos(cos(x))
#'
#' @param x row matrix
#'
#'  @return n x n similarity matrix
#'  @export
#'
#'  @examples
cosine <- function(x){

  s <- tcrossprod(x) / (sqrt(tcrossprod(rowSums(x^2))))

  diag(s) <- 1

  s

}
