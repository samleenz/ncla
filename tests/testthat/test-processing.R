library(ncla)
# make toy data
ensg <- readRDS(testthat::test_path("ensembl_genes.rds"))
print(getwd())

df <- matrix(data = pmax(rnorm(1000*25, 100, 50), 0), ncol = 10000, nrow = 25, dimnames = list(LETTERS[1:25], sample(ensg, 10000)))


my_geneset <- sample(colnames(df), 1500)

test_that("custom gene-sets work", {
  expect_error(subsetGenes(df, "Custom__missing_geneset"))
  expect_equal(ncol(subsetGenes(df, "Custom__my_geneset")), 1500)
})

test_that("geneset name is valid", {
  expect_(subsetGenes(df, "Hallmarks")),
  expect_error(subsetGenes(df, "fail")),
  expect_equal(subsetGenes(df, "Full"),df),
  expect_error(subsetGenes(df, "haLlMarKs"))
})
