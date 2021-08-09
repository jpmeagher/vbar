test_that("individual specific latent trait variance", {
  D <- 15
  L <- 4
  W <- array(rnorm(D*L), c(D, L))
  W_outer <- apply(W, 1, function(x) x %*% t(x)) %>%
    array(c(L, L, D))
  lambda <- rgamma(D, 1, 1)
  tau <- runif(L, 0, 0.5)

  expect_equal(
    W_outer[, , 1], W[1, ] %*% t(W[1, ])
  )

  tmp <- diag(1 / tau^2)
  for (i in 1:D) {
    tmp <- tmp + lambda[i] * W_outer[, , i]
  }

  expect_equal(
    individual_specific_latent_trait_precision(
      auxiliary_trait_precision_vector = lambda,
      loading_outer_product_expectation = W_outer,
      within_taxon_standard_deviation = tau,
      perform_checks = TRUE
    ),
    tmp
  )
})
