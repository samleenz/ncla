#' add nearest neighbours
#'
#' Add the nearest reference sample to the query data in a SummarizedExperiment
#' object.
#'
#' @param se SummarisedExperiment object
#' @param res the output of [ncla()], the samples in the query column should
#'   have identical names as in the `se` object
#' @param colname character (default `nearest_reference`), new colname to add to
#'   the `colData()` slot
#'
#' @return a summarised experiment object with the nearest reference sample for
#'   each query added to the colData slot
#' @export
#'
#' @examples
add_NN <- function(se, res, colname = "nearest_reference"){
  res_filt <- res[res$rank == 1, c("query", "reference")]

  # make a dict with
  #     key -> query
  #     val -> reference (rank 1)
  ref_map <- res_filt$reference
  names(ref_map) <- res_filt$query

  if(colname %in% colnames(SummarizedExperiment::colData(se))){
    stop(glue::glue("{colname} already in colData(se), choose a new name"))
  }

  SummarizedExperiment::colData(se)[, colname] <- ref_map[colnames(se)] |>
    unname()

  se
}
