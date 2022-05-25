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
    d <- philentropy::distance(x, "manhattan", mute.message = TRUE, use.row.names = TRUE)
  } else if(f == "L2"){
    d <- philentropy::distance(x, "euclidean", mute.message = TRUE, use.row.names = TRUE)
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
      mute.message = TRUE,
      use.row.names = TRUE
      )
  } else if(f == "bhattacharyya"){
    d <- philentropy::distance(
      x,
      "bhattacharyya",
      est.prob = "empirical",
      mute.message = TRUE,
      use.row.names = TRUE
      )
  } else if(f == "F5"){ # fractional minkowksi distnace, p = 0.5
    d <- philentropy::distance(
      x,
      "minkowski",
      p = 0.5,
      mute.message = TRUE,
      use.row.names = TRUE
      )
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




#' (internal) cross distance
#'
#' @param x length p vector
#' @param y matrix with p cols
#' @param f
#' @param g
#'
#' @return
#'
#' @examples
.xdissim <- function(x, y, f, g = NULL){
  # for each row in y, calculate f() with x
  purrr::map(
    rownames(y),
    ~ rbind(y[.x, ], x)
  ) |>
    purrr::map(dissim, f = f, g = g) |>
    purrr::map_dbl(2) # grab the top-right of the 2x2 distance
  #                     matrix that is returned
}

#' cross distance
#'
#' Given two matrices, x and y, find the distance between each row of x to each row of y
#'
#' @param x a row matrix (or coercable dataframe)
#' @param y a row matrix (or coercable dataframe)
#' @param f function to use, see dissim
#' @param g gene space for tsd, see dissim
#'
#' @return
#' @export
#'
#' @examples
xdissim <- function(x, y, f, g = NULL){
  # given matrix x, get dissim of each row of x to each row of y
  # return as nrow(x) x nrow(y) matrix
  raw_d <- furrr::future_map(
    rownames(x),
    ~ .xdissim(x[.x,], y, f, g)
  )

  # convert list of distances to matrix
  matrix(unlist(raw_d), ncol = nrow(y), byrow = TRUE, dimnames = list(rownames(x), rownames(y)))
}

#' k nearest samples
#'
#' Given a matrix of cross distances, return the k nearest neighbours
#' for each of the queries (rows)
#'
#' @param xd
#' @param k
#'
#' @return a `nrow(xd)` length list with k named + sorted elements per list
#' @export
#'
#' @examples
kns <- function(xd, k = 1){
  # take cross dist matrix and return k nearest samples
  neighbs <- purrr::map(
    rownames(xd),
    ~ rank(xd[.x, ], ties.method = "min")
  ) |>
    purrr::map(~ .x[.x <= k]) |> # keep only those below our neighbour cutoff
    purrr::map(sort)# |>

  names(neighbs) <- rownames(xd)

  neighbs
}

#' kns tidy
#'
#' transform ks output into a tidy table
#'
#' @param x output of kns()
#'
#' @return 3 column tibble (query, reference, rank)
#' @export
#'
#' @examples
kns_tidy <- function(x){
  x |>
    purrr::map_df(tibble::enframe, name = "reference", value = "rank", .id = "query") |>
    dplyr::select(query, reference, rank)
}
