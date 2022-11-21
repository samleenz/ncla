#' Batch correct
#'
#' Wrapper for batch correction of two datasets
#'
#' Make sure you have the package needed install (limma or batchelor)
#'
#' @param x datasset one: sample x gene matrix
#' @param y datasset one: sample x gene matrix
#' @param method oneof c("MNN", "limma")
#' @param prop.k see batchelor::fastMNN
#'
#' @return legnth two list of batch corrected datasets
#' @export
#'
#' @examples
batchCorrect <- function(x, y, method = "limma", prop.k = 0.05){
  if (!requireNamespace("pkg", quietly = TRUE)) {
    stop(
      "Package \"pkg\" must be installed to use this function.",
      call. = FALSE
    )
  }


  # check x and y are dataframes

  # make feature space identical
  features <- intersect(colnames(x), colnames(y))
  x <- x[, features]
  y <- y[, features]


  # prevent duplicate sample names, store original for later
  rn_x <- rownames(x)
  rownames(x) <- paste("one", rownames(x), sep = "_")
  rn_y <- rownames(y)
  rownames(y) <- paste("two", rownames(y), sep = "_")

  if(method == "limma"){

    # make combined data
    dat <- rbind(x, y) %>%
      t()

    dat_bc <- transformData(dat, "log2") %>%
      limma::removeBatchEffect(
        batch = c(rep("one", length(rn_x)), rep("two", length(rn_y)))
        )

    x_bc <- .getBatchLimma(dat_bc, "one")
    y_bc <- .getBatchLimma(dat_bc, "two")

  } else if(method == "MNN"){
    # print("did MNN")
    mnn <- batchelor::fastMNN(
      t(transformData(x, "log2")),
      t(transformData(y, "log2")),
      prop.k = prop.k
    )

    x_bc <- .getBatchMNN(mnn, "one")
    y_bc <- .getBatchMNN(mnn, "two")

  } else{
    stop(paste(method, "is not a defined batch correction method"))
  }


  return(
    list(
      x = x_bc,
      y = y_bc
    )
  )
}


#' get batch
#'
#' internal function to extract individual batch from limma corrected df
#'
#' @param df
#' @param batch
#'
#' @return
#'
#' @examples
.getBatchLimma <- function(df, batch){
  if(! any(startsWith(colnames(df), batch))){
    stop("Batch is not present in df")
  }

  x <- df[, startsWith(colnames(df), batch)]
  colnames(x) <- stringr::str_remove(colnames(x), paste0(batch, "_"))
  t(x)
}

#' get batch
#'
#' internal function to extract individual batch from MNN corrected df
#'
#' @param sce
#' @param batch
#'
#' @return
#'
#' @examples
.getBatchMNN <- function(sce, batch){
  # MNN uses numerical numbering so quick lookup table to smooth over
  batches <- c(one = 1, two = 2)

  df <- sce[, sce$batch == batches[[batch]]]

  x <- df %>%
    SummarizedExperiment::assay("reconstructed") %>%
    as.data.frame()

  colnames(x) <- stringr::str_remove(colnames(x), paste0(batch, "_"))
  t(x)

}
