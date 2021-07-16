test_that("values simulated for OU process", {
  tr <- scale_phylo(ape::rtree(100))

  z <- simulate_phylogenetic_ou(tr)
  checkmate::expect_numeric(
    z, any.missing = FALSE, len = 100, names = "named"
    )
  expect_equal(
    names(z), tr$tip.label
  )

  z <- simulate_phylogenetic_ou(
    tr, internal = TRUE
    )
  checkmate::expect_numeric(
    z, any.missing = FALSE, len = 200-1, names = "named"
  )
  expect_equal(
    names(z[1:100]), tr$tip.label
  )
  expect_equal(
    names(z[-(1:100)]), paste("Node", 1:99, sep = "")
  )
})

