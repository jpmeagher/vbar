test_that("latent_trait outer product", {
  L <- 5
  z <- rnorm(L)
  S <- stats::rWishart(1, df = L, Sigma = diag(L)/ L^2)[,,1]

  expect_equal(
    latent_trait_outer_product_expectation(z_n_tilde = z, S_n_z = S),
    S + z %*% t(z)
  )

})
