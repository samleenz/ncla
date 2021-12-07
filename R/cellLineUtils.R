#' clean cell IDs
#'
#' given a character vector of cell-line IDs, return a standardised version.
#'
#' all in UPPER case, no punctuation, and some other idiosyncrasies removed
#'
#' @param x character vector of cell IDs
#' @param dataset character denoting which cell line the IDs `x` are from
#'     One of
#'       - nci60Full
#'       - klijn2015Full
#'       - ccleFull
#'
#' @return
#' @export
#'
#' @examples
cleanCellID <- function(x, dataset){
  # check source is a valid column name
  if(! dataset %in% functions$clean_IDs$dataset){
    stop(paste(dataset, "is not a valid data-set selection"))
  }

  df <- functions$clean_IDs[functions$clean_IDs$dataset == dataset, ]
  df <- df[match(x, df$original_id), ]

  df[["cellid"]]
}
