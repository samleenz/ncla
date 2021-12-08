#' transformData: apply a transformtation to row matrix
#'
#' @param x the row matrix
#' @param f transformation to apply, one of:
#'   - none: return the same matrix
#'   - log2: log-2 transform with pseudo-count of one
#'   - scale: scale and center
#'   - pcaN  : PCA where N PCs are kept
#'
#'
#' @return the transformed data matrix
#' @export
#'
#' @examples
transformData <- function(x, f){
  if(f == "none"){
    tr <- x
  } else if(f == "log2"){
    tr <- log2(x  + 1)
  } else if(f == "scale"){
    # scale each sample
    tr <- t(scale(t(x)))
  } else if(startsWith(f, "pca")){
    nTop <- stringr::str_remove(f, "pca")
    if(is.na(as.numeric(nTop))){
      stop("culmulative var % for PCA cannot be coerced to a numeric")
    } else if(as.numeric(nTop) > min(dim(x))){
      stop("Too many PCs requested, cannot be more than shortest dim of x")
    } else {
      tr <- pcaTopN(x, as.numeric(nTop))
    }
  } else if(f == "umap"){
    if(nrow(x) < umap::umap.defaults$n_neighbors){
      umap_settings <- umap::umap.defaults
      umap_settings$n_neighbors <- nrow(x) - 1
    } else{
      umap_settings <- umap::umap.defaults
    }
    tr <- umap::umap(x, config = umap_settings)[["layout"]]
  } else {
    stop("f is not a valid option")
  }

  tr
}

pcaTopN <- function(x, i){
  pca <- prcomp(x, scale. = TRUE)

  PCs <- pca$x[, 1:i]

  return(PCs)
}

#' Subset genes
#'
#' subset a gene-expression matrix with a defined geneset
#'
#' @param x row matrix where columns are ENSG labelled gens
#' @param f character name of geneset to use
#'   - Halmarks
#'   - L1000
#'   - Stable
#'
#' @return
#' @export
#'
#' @examples
subsetGenes <- function(x, f){
  if(f == "Hallmarks"){
    genes <- gs_hallmark[gs_hallmark %in% colnames(x)]
  } else if(f == "L1000"){
    genes <- gs_l1000[gs_l1000 %in% colnames(x)]
  } else if(f == "Stable"){
    genes <- gs_l1000[gs_l1000 %in% colnames(x)]
  } else if(startsWith(f, "Custom")) {
    # f is the form 'Custom__<name>' where name should point to a vector of
    # gene names in the environment
    nme <- stringr::str_split(f, "__", simplify = T)[2]
    g <- mget(nme, ifnotfound = FALSE, envir = .GlobalEnv)[[1]]
    if(isFALSE(g)){
      # if nme not found in namespace
      message(paste(nme,  "did not match, returning x"))
      genes <- colnames(x)
    } else {
      genes <- g[g %in% colnames(x)]
    }
  } else {
    message(paste(f,  "did not match, returning x"))
    genes <- colnames(x)

    apply
  }

  x[, genes]
}
