test_that("individual specific latent trait precision", {
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
    compute_individual_specific_latent_trait_precision(
      auxiliary_trait_precision_vector = lambda,
      loading_outer_product_expectation = W_outer,
      within_taxon_standard_deviation = tau,
      perform_checks = TRUE
    ),
    tmp
  )
})

test_that("individual specific latent trait expectation", {
  D <- 15
  L <- 4
  W <- array(rnorm(D*L), c(D, L))
  W_outer <- apply(W, 1, function(x) x %*% t(x)) %>%
    array(c(L, L, D))
  lambda <- rgamma(D, 1, 1)
  tau <- runif(L, 0, 0.5)

  inv_S_z <- compute_individual_specific_latent_trait_precision(
    auxiliary_trait_precision_vector = lambda,
    loading_outer_product_expectation = W_outer,
    within_taxon_standard_deviation = tau,
    perform_checks = TRUE
  )

  z <- rnorm(L)
  x <- W %*% z + rnorm(D, sd = sqrt(1 / lambda))

  f <- rnorm(L)

  S_z <- solve(inv_S_z)
  tmp <- S_z %*% (
    (t(W) %*% diag(lambda) %*% x) + (diag(1 / (tau^2)) %*% f)
  )

  expect_equal(
    compute_individual_specific_latent_trait_expectation(
      auxiliary_trait = x,
      loading = W,
      taxon_specific_latent_trait = f,
      auxiliary_trait_precision_vector = lambda,
      individual_specific_latent_trait_precision = inv_S_z,
      within_taxon_standard_deviation = tau,
      perform_checks = TRUE
    ),
    c(tmp)
  )
})

test_that("terminal taxon-specific latent trait precision", {
  N <- 10
  L <- 4
  tau <- runif(L, 0, 0.5)
  eta <- runif(L, 0, 0.5)

  tmp <- N * (1 / tau^2) + (1 / eta^2)

  expect_equal(
    compute_terminal_taxon_specific_latent_trait_precision(
      N = N,
      within_taxon_standard_deviation = tau,
      conditional_standard_deviation = eta
    ),
    tmp
  )
})


