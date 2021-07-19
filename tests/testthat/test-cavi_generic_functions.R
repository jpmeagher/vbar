test_that("Gaussian outer product", {
  L <- 5
  z <- rnorm(L)
  S <- stats::rWishart(1, df = L, Sigma = diag(L)/ L^2)[,,1]

  expect_equal(
    gaussian_outer_product_expectation(expected_value = z, covariance_matrix = S),
    S + z %*% t(z)
  )

})

