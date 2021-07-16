test_that("phylogenetic ou hyperparameters are computed correctly", {
  n <- 100
  alpha <- 10
  ell <- 1/2

  tr <- ape::rtree(n)
  hyper <- phylogenetic_ou_hyperparameters(
    tr,
    amplitude = alpha, length_scale = ell
    )
  checkmate::expect_matrix(
    hyper, any.missing = FALSE, nrows = 2*n - 1, ncols = 2, col.names = "named"
  )

  D <- ape::dist.nodes(tr)
  K <- ou_kernel(d = D, amplitude = alpha, length_scale = ell)
  tmp_hyper <- matrix(nrow = 2*n - 1, ncol = 2)
  tmp_hyper[n+1, ] <- c(1, K[n+1, n+1])
  for(i in seq_along(tr$edge.length)) {
    tmp_hyper[tr$edge[i, 2], 1] <- K[tr$edge[i, 1], tr$edge[i, 2]] / K[tr$edge[i, 1], tr$edge[i, 1]]
    tmp_hyper[tr$edge[i, 2], 2] <- K[tr$edge[i, 2], tr$edge[i, 2]] - tmp_hyper[tr$edge[i, 2], 1] * K[tr$edge[i, 2], tr$edge[i, 1]]
  }
  colnames(tmp_hyper) <- c("parent_weight", "variance")

  expect_equal(
    hyper, tmp_hyper
  )
})
