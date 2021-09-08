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
      precision_vector = lambda,
      loading_outer_product_expectation = W_outer,
      within_taxon_amplitude = tau,
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
    precision_vector = lambda,
    loading_outer_product_expectation = W_outer,
    within_taxon_amplitude = tau,
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
      precision_vector = lambda,
      individual_specific_latent_trait_precision = inv_S_z,
      within_taxon_amplitude = tau,
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
      within_taxon_amplitude = tau,
      conditional_standard_deviation = eta
    ),
    tmp
  )
})

test_that("terminal taxon-specific latent trait expectation", {
  N <- 111
  L <- 4
  tau <- runif(L, 0, 0.5)
  eta <- runif(L, 0, 0.5)
  nu <- runif(L, 0.5, 1)

  f_pa <- rnorm(L)
  f <- nu * f_pa + rnorm(L, sd = eta)
  Z <- sweep(
    matrix(rnorm(N*L, sd = tau), nrow = N, ncol = L, byrow = TRUE),
    2, f, "+"
    )

  inv_S_f <- compute_terminal_taxon_specific_latent_trait_precision(
    N = N,
    within_taxon_amplitude = tau,
    conditional_standard_deviation = eta
  )
  S_f <- 1 / inv_S_f

  tmp <- S_f * (colSums(Z %*% diag(1 / tau^2)) + (nu * f_pa / (eta^2)))

  expect_equal(
    compute_terminal_taxon_specific_latent_trait_expectation(
      individual_specific_latent_traits = Z,
      within_taxon_amplitude = tau,
      parent_taxon_latent_trait = f_pa,
      conditional_expectation_weight = nu,
      conditional_standard_deviation = eta,
      latent_trait_precision = inv_S_f,
      perform_checks = TRUE
    ),
    tmp
  )
})

test_that("internal taxon-specific latent trait precision", {
  N <- 2
  L <- 4

  eta <- runif(L, 0, 0.5)
  eta_ch <- matrix(runif(N * L, 0, 0.5), nrow = N, ncol = L)
  nu_ch <- matrix(runif(N * L, 0.5, 1), nrow = N, ncol = L)

  tmp <- (nu_ch / eta_ch) %>%
    magrittr::raise_to_power(2) %>%
    apply(2, sum) %>%
    magrittr::add(1 / eta^2)

  expect_equal(
    compute_internal_taxon_specific_latent_trait_precision(
      child_taxa_conditional_expectation_weights = nu_ch,
      child_taxa_conditional_standard_deviations = eta_ch,
      conditional_standard_deviation = eta
    ),
    tmp
  )
})

test_that("internal taxon-specific latent trait expectation", {
  N <- 2
  L <- 4

  eta <- runif(L, 0, 0.5)
  nu <- runif(L, 0.5, 1)
  eta_ch <- matrix(runif(N * L, 0, 0.5), nrow = N, ncol = L)
  nu_ch <- matrix(runif(N * L, 0.5, 1), nrow = N, ncol = L)

  f_pa <- rnorm(L)
  f <- nu * f_pa + rnorm(L, sd = eta)
  f_ch <- sweep(
    matrix(rnorm(N*L, sd = nu_ch), nrow = N, ncol = L, byrow = TRUE),
    2, f, "+"
  )

  inv_S_f <- compute_internal_taxon_specific_latent_trait_precision(
    child_taxa_conditional_expectation_weights = nu_ch,
    child_taxa_conditional_standard_deviations = eta_ch,
    conditional_standard_deviation = eta
  )

  tmp <- (nu_ch * f_ch / (eta_ch^2)) %>%
    apply(2, sum) %>%
    magrittr::add(
      nu * f_pa / (eta^2)
    ) %>%
    magrittr::divide_by(
      inv_S_f
    )

  expect_equal(
    compute_internal_taxon_specific_latent_trait_expectation(
      child_taxa_latent_traits = f_ch,
      child_taxa_conditional_expectation_weights = nu_ch,
      child_taxa_conditional_standard_deviations = eta_ch,
      parent_taxon_latent_trait = f_pa,
      conditional_expectation_weight = nu,
      conditional_standard_deviation = eta,
      latent_trait_precision = inv_S_f,
      perform_checks = TRUE
    ),
    tmp
  )
})
