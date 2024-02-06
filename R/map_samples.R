#' map_samples
#'
#' For a given set of `query` samples, perform dissimilarity mapping against the
#' supplied `reference` using the dissimilarity measure specified by `f`. This
#' function enforces that the features are the same, e.g. `query = x[, genes]`
#' and `reference = y[, genes]`.
#'
#' Internally, this function uses foreach with the
#' ability to run in parallel, see [doParallel::registerDoParallel()] for how to
#' register a parallel backend.
#'
#' @param query row matrix with samples to query
#' @param reference row matrix with samples to use as the reference
#' @param f The dissimilarity function to use, see [dissim()]
#'
#' @return a data.frame
#' @export
#'
#' @examples
map_samples <- function(query, reference, f){
  if(! all(colnames(query) == colnames(reference))){
    stop("Columns of query and reference must be identical")
  }

  foreach(cl = rownames(query), .combine = "rbind", .errorhandling = "remove") %dopar% {
    xyMat <- rbind(xi = query[cl, ], reference) %>%
      as.matrix()

    di <- sapply(2:nrow(xyMat), function(i){
      ncla::dissim(
        xyMat[c(1, i), ],
        f
      )[1,2] # [1,2] is to select the top-right of the 2x2 distance matrix
    })

    data.frame(
      query = cl,
      reference = rownames(xyMat)[-1],
      dist = di,
      rank = rank(di, ties.method = "min")
    )
  }
}
