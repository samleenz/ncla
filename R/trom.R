#' calcAssociatedGenes
#'
#' given a data-set, z-score transform and get associated genes as per TROM method
#'
#' @param x a n x m data frame, normalised gene expression with *n* samples in rows
#' @param z z-score threshold for associated genes
#'
#' @return length(n) list with associated genes for each sample
#' @export
#'
#' @examples
calcAssociatedGenes <- function(x, z = 1.5){
  x_z <- as.matrix(scale(x, center = TRUE, scale = TRUE))

  ass_genes <- apply(
    x_z,
    1,
    function(x){
      x1 <- x[x > z]
      names(x1)
    },
    simplify = F
  )

  ass_genes
}

trom2 <- function(x, z = 1.5){

  # get associated genes for each sample
  ass_genes <- calcAssociatedGenes(x, z)
  # 🍑🧬 is a list length(ncol(samples))

  # create square matrix for output
  s <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  dimnames(s) <- list(
    names(ass_genes),
    names(ass_genes)
  )

  # calculate overlap p-values
  num_tests <- sum(upper.tri(s)) + length(diag(s)) # this is the number of tests we do
  for(i in seq(nrow(x))){
    for(j in seq(nrow(x))){
      if(i <= j){
        # num_tests <- num_tests + 1

        # Run Fisher's eaxct test to test enrichment of overlap between
        # ass_genes[[i]] and ass_genes[[j]]
        # https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
        g1  <-  length(ass_genes[[i]])
        g2  <-  length(ass_genes[[j]])
        overlap <-  length(intersect(ass_genes[[i]], ass_genes[[j]]))
        total   <- ncol(x)

        fisher_p <- fisher.test(
          matrix(
            c(
              overlap,
              g2 - overlap,
              g1 - overlap,
              total - g2 - g1 + overlap
            ),
            nrow = 2,
            ncol = 2
          ),
          alternative='greater'
        )$p.value

        # adjust p-values by number of tests before assigning
        s[i, j] <- p.adjust(fisher_p, "BH", n = num_tests)
        s[j, i] <- p.adjust(fisher_p, "BH", n = num_tests)
      }
    }
  }

  # turn p-values into TROM scores
  # -log10(adjusted_pvalue)
  # perfect matches (p-value = 0 get 300)
  s_t <- - log10(s)
  s_t[is.infinite(s_t)] <- 300

  return(s_t)
}


#' TROM from associated genes
#'
#' given two sets of associated genes, calculate the TROM score
#'
#' @param ag1 first list of associated genes
#' @param ag2 second list of associated genes
#' @param total total number of genes in the dataset
#' @param n_test total number of tests for multiple testing correction
#'
#' @return
#' @export
#'
#' @examples
trom_ag <- function(ag1, ag2, total, n_test = 1){
  g1 <- length(ag1)
  g2 <- length(ag2)
  overlap <- length(intersect(ag1, ag2))
  total <- total # just for readiblity sake

  fisher_p <- fisher.test(
    matrix(
      c(
        overlap,
        g2 - overlap,
        g1 - overlap,
        total - g2 - g1 + overlap
      ),
      nrow = 2,
      ncol = 2
    ),
    alternative='greater'
  )$p.value

  -log10(p.adjust(fisher_p, "BH", n = n_test))

}
