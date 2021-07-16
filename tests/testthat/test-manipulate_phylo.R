test_that("phylogeny is scaled correctly", {
  tr <- ape::rtree(100)
  tr1 <- scale_phylo(tr)
  phytools::nodeHeights(tr)
  expect_equal(
    max(phytools::nodeHeights(tr1)), 1
  )
})
