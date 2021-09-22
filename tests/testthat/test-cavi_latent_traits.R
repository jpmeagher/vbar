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

test_that("individual specific latent trait contribution to the elbo", {
  P <- 4
  L <- 4
  tn <- c("ord", "nom", "con", "fvt")
  tt <- factor(tn, levels = c("ord", "nom", "con", "fvt"))
  ind_mt <- list(
    ord = 2L, nom = 3L, con = 4L, fvt = 5:36
  )
  mt <- cbind(synthetic_traits[, 1:4], fvt = t(simplify2array(synthetic_traits$fvt)))
  N <- nrow(mt)
  K <- c(4, 3, NA, NA)
  ind_at <- list(
    ord = 1L, nom = 1 + 1:3, con = 5L, fvt = 5 + 1:32
  )
  D_p <- sum(sapply(ind_at, length))

  cat <- list(
    ord = NA,
    nom = factor(levels(mt$nom)),
    con = NA, fvt = NA
  )
  gamma <- list(
    ord = c(-Inf, 0, 1, 2, Inf),
    nom = NA, con = NA, fvt = NA
  )

  g <- list(
    ord = ordinal_link, nom = nominal_link, con = function(x) x, fvt = function(x) exp(x)
  )
  g_inv <- list(
    ord = function(y){
      ordinal_inverse_link(
        y, cut_off_points = gamma$ord,
        mu = rep(0, N), return_expectation = FALSE
      )
    },
    nom = function(y){
      nominal_inverse_link(
        y, mu = matrix(0, N, length(cat$nom)),
        n_samples = 1000, return_expectation = FALSE
      )
    },
    con = function(y) y,
    fvt = function(y) log(data.matrix(y))
  )

  meta <- specify_manifest_trait_metadata(
    n_traits = P, trait_names = tn, trait_type = tt,
    trait_levels = K,
    manifest_trait_index = ind_mt, auxiliary_trait_index = ind_at,
    link_functions = g,
    inverse_link_functions = g_inv,
    cut_off_points = gamma, categories = cat,
    manifest_trait_df = mt,
    perform_checks = TRUE
  )

  C_w <- diag(D_p)
  x <- seq(0, 1, length.out = length(ind_at[[4]]))
  d <- abs(outer(x, x, "-"))
  ell <- 1 / (2 * pi)
  C_w[ind_at[[4]], ind_at[[4]]] <-  (exp_quad_kernel(d, 1, ell) + (1e-6 * diag(length(ind_at[[4]])))) / (1 + 1e-6)

  ph <- vbar::synthetic_trait_model_specification$phylogeny
  S <- length(ph$tip.label)

  plvm <- initialise_plvm(
    manifest_trait_df = mt, metadata = meta, phy = ph,
    L = L,
    loading_prior_correlation = C_w,
    auxiliary_traits = NULL,
    precision_prior_shape = 1, precision_prior_rate = 0.01,
    precision = NULL,
    ard_precision = NULL,
    ard_shape = 1, ard_rate = 1,
    loading = NULL, method = "random",
    within_taxon_amplitude = NULL,
    heritable_amplitude = NULL,
    length_scale = 2,
    perform_checks = TRUE
  )

  npt <- table(
    mt$taxon_id
  )[ph$tip.label] %>% as.numeric()

  elbo <- compute_individual_specific_latent_trait_elbo(
    individual_specific_latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
    taxon_id = mt$taxon_id, phy = ph,
    terminal_taxon_specific_latent_trait_expectation = plvm$taxon_specific_latent_trait_expectation[1:S,],
    individual_specific_latent_trait_covariance = plvm$individual_specific_latent_trait_covariance,
    individual_specific_latent_trait_outer_product_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
    terminal_taxon_latent_trait_outer_product_expectation = plvm$taxon_specific_latent_trait_outer_product_expectation[, , 1:S],
    within_taxon_amplitude = plvm$within_taxon_amplitude,
    perform_checks = TRUE
  )

  tmp1 <- apply(plvm$individual_specific_latent_trait_outer_product_expectation, c(1, 2), sum)
  tmp2 <- sweep(
    plvm$taxon_specific_latent_trait_outer_product_expectation[, , 1:S],
    3, npt, "*"
  ) %>%
    apply(c(1, 2), sum)

  expect_equal(
    elbo,
    - 0.5 * N * determinant(diag(plvm$within_taxon_amplitude^2))$modulus[1] -
      (0.5 * sum(diag((tmp1 + tmp2) %*% diag(plvm$within_taxon_amplitude^-2)))) +
      sum(sapply(
        1:S, function(i){
          t(plvm$taxon_specific_latent_trait_expectation[i, ]) %*%
            diag(plvm$within_taxon_amplitude^-2) %*%
            colSums(
              plvm$individual_specific_latent_trait_expectation[mt$taxon_id == ph$tip.label[i], ]
            )
        })) +
      0.5 * N * determinant(plvm$individual_specific_latent_trait_covariance)$modulus[1]
  )
})

test_that("taxon specific latent trait contribution to the elbo", {
  P <- 4
  L <- 4
  tn <- c("ord", "nom", "con", "fvt")
  tt <- factor(tn, levels = c("ord", "nom", "con", "fvt"))
  ind_mt <- list(
    ord = 2L, nom = 3L, con = 4L, fvt = 5:36
  )
  mt <- cbind(synthetic_traits[, 1:4], fvt = t(simplify2array(synthetic_traits$fvt)))
  N <- nrow(mt)
  K <- c(4, 3, NA, NA)
  ind_at <- list(
    ord = 1L, nom = 1 + 1:3, con = 5L, fvt = 5 + 1:32
  )
  D_p <- sum(sapply(ind_at, length))

  cat <- list(
    ord = NA,
    nom = factor(levels(mt$nom)),
    con = NA, fvt = NA
  )
  gamma <- list(
    ord = c(-Inf, 0, 1, 2, Inf),
    nom = NA, con = NA, fvt = NA
  )

  g <- list(
    ord = ordinal_link, nom = nominal_link, con = function(x) x, fvt = function(x) exp(x)
  )
  g_inv <- list(
    ord = function(y){
      ordinal_inverse_link(
        y, cut_off_points = gamma$ord,
        mu = rep(0, N), return_expectation = FALSE
      )
    },
    nom = function(y){
      nominal_inverse_link(
        y, mu = matrix(0, N, length(cat$nom)),
        n_samples = 1000, return_expectation = FALSE
      )
    },
    con = function(y) y,
    fvt = function(y) log(data.matrix(y))
  )

  meta <- specify_manifest_trait_metadata(
    n_traits = P, trait_names = tn, trait_type = tt,
    trait_levels = K,
    manifest_trait_index = ind_mt, auxiliary_trait_index = ind_at,
    link_functions = g,
    inverse_link_functions = g_inv,
    cut_off_points = gamma, categories = cat,
    manifest_trait_df = mt,
    perform_checks = TRUE
  )

  C_w <- diag(D_p)
  x <- seq(0, 1, length.out = length(ind_at[[4]]))
  d <- abs(outer(x, x, "-"))
  ell <- 1 / (2 * pi)
  C_w[ind_at[[4]], ind_at[[4]]] <-  (exp_quad_kernel(d, 1, ell) + (1e-6 * diag(length(ind_at[[4]])))) / (1 + 1e-6)

  ph <- vbar::synthetic_trait_model_specification$phylogeny
  S <- length(ph$tip.label)

  plvm <- initialise_plvm(
    manifest_trait_df = mt, metadata = meta, phy = ph,
    L = L,
    loading_prior_correlation = C_w,
    auxiliary_traits = NULL,
    precision_prior_shape = 1, precision_prior_rate = 0.01,
    precision = NULL,
    ard_precision = NULL,
    ard_shape = 1, ard_rate = 1,
    loading = NULL, method = "random",
    within_taxon_amplitude = NULL,
    heritable_amplitude = NULL,
    length_scale = 2,
    perform_checks = TRUE
  )

  to <- ph$edge[ape::postorder(ph), ]

  elbo <- compute_taxon_specific_latent_trait_elbo(
    taxon_specific_latent_trait_expectation = plvm$taxon_specific_latent_trait_expectation,
    taxon_specific_latent_trait_outer_product_expectation = plvm$taxon_specific_latent_trait_outer_product_expectation,
    taxon_specific_latent_trait_covariance = plvm$taxon_specific_latent_trait_covariance,
    phy = ph,
    phylogenetic_gp = plvm$phylogenetic_GP,
    perform_checks = TRUE
  )

  tmp1 <- sum(sapply(
    1:(2 * S - 2), function(i){
      sum(diag((plvm$taxon_specific_latent_trait_outer_product_expectation[, , to[i, 2]] +
              (diag(plvm$phylogenetic_GP[to[i, 2], "weight", ]) %*%
              plvm$taxon_specific_latent_trait_outer_product_expectation[, , to[i, 1]] %*%
              diag(plvm$phylogenetic_GP[to[i, 2], "weight", ]))
              ) %*% diag(plvm$phylogenetic_GP[to[i, 2], "sd", ]^-2)))
    }
  ))
  tmp2 <- sum(diag(plvm$taxon_specific_latent_trait_outer_product_expectation[, , S + 1] %*%
                     diag(plvm$phylogenetic_GP[S+1, "sd", ]^-2)))

  expect_equal(
    elbo,
    - 0.5 * sum(sapply(
      1:(2 * S - 1), function(i){
        determinant(diag(plvm$phylogenetic_GP[i, "sd", ]^2))$modulus[1]
      }
    )) -
      0.5 * (tmp1 + tmp2) +
      sum(sapply(
        1:(2 * S - 2), function(i){
          t(plvm$taxon_specific_latent_trait_expectation[to[i, 1], ]) %*%
            diag(plvm$phylogenetic_GP[to[i, 2], "weight", ]) %*%
            diag(plvm$phylogenetic_GP[to[i, 2], "sd", ]^-2) %*%
            plvm$taxon_specific_latent_trait_expectation[to[i, 2], ]
        }
      )) +
      0.5 * sum(sapply(
        1:(2 * S - 1), function(i){
          determinant(exp(1) * diag(plvm$taxon_specific_latent_trait_covariance[i, ]))$modulus[1]
        }
      ))
  )
})

