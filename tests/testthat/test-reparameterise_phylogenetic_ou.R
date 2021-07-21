test_that("phylogenetic ou hyperparameters are computed correctly", {
  S <- 10
  alpha <- 2
  ell <- 3
  sigma <- 0.5

  tr <- ape::rtree(S)
  par <- reparameterise_phylogenetic_ou(
    tr,
    heritable_amplitude = alpha, length_scale = ell,
    environmental_amplitude = sigma
    )
  checkmate::expect_matrix(
    par, any.missing = FALSE, nrows = 2*S - 1, ncols = 2, col.names = "named"
  )

  D <- ape::dist.nodes(tr)
  K <- ou_kernel(d = D, amplitude = alpha, length_scale = ell)
  for (i in 1:S) K[i, i] <-  K[i, i] + sigma^2
  tmp_par <- matrix(nrow = 2*S - 1, ncol = 2)
  tmp_par[S+1, ] <- c(1, K[S+1, S+1])
  for(i in seq_along(tr$edge.length)) {
    tmp_par[tr$edge[i, 2], 1] <- K[tr$edge[i, 1], tr$edge[i, 2]] / K[tr$edge[i, 1], tr$edge[i, 1]]
    tmp_par[tr$edge[i, 2], 2] <- K[tr$edge[i, 2], tr$edge[i, 2]] - tmp_par[tr$edge[i, 2], 1] * K[tr$edge[i, 2], tr$edge[i, 1]]
  }
  tmp_par[, 2] <- sqrt(tmp_par[, 2])
  colnames(tmp_par) <- c("weight", "sd")

  expect_equal(
    par[, 1], tmp_par[, 1]
  )
  expect_equal(
    par[2], tmp_par[2]
  )
})
