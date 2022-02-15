test_that("Loading initialisation works", {
  P <- 4
  L <- 4
  D_p <-
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

  X <- initialise_auxiliary_traits(
    n_traits = nrow(meta),
    manifest_trait_df = mt,
    trait_names = meta$trait_names,
    trait_type = meta$trait_type,
    trait_levels = meta$trait_levels,
    manifest_trait_index = meta$manifest_trait_index,
    auxiliary_trait_index = meta$auxiliary_trait_index,
    inverse_link_functions = meta$inverse_link_functions,
    auxiliary_traits = NULL,
    perform_checks = TRUE
  )

  C_w <- diag(D_p)
  x <- seq(0, 1, length.out = length(ind_at[[4]]))
  d <- abs(outer(x, x, "-"))
  ell <- 1 / (2 * pi)
  C_w[ind_at[[4]], ind_at[[4]]] <-  (exp_quad_kernel(d, 1, ell) + (1e-6 * diag(length(ind_at[[4]])))) / (1 + 1e-6)

  alpha <- rgamma(L, shape = 1, rate = 1)

  W <- initialise_loading(
    D_prime = D_p, L = L,
    ard_precision = alpha, loading_prior_correlation = C_w,
    loading = NULL, method = "random",
    auxiliary_traits = NULL,
    perform_checks = TRUE
  )

  checkmate::expect_matrix(
    W, mode = "numeric", any.missing = FALSE,
    nrows = D_p, ncols = L
  )

  W <- initialise_loading(
    D_prime = D_p, L = L,
    ard_precision = alpha, loading_prior_correlation = C_w,
    loading = NULL, method = "pca",
    auxiliary_traits = X,
    perform_checks = TRUE
  )

  pca <- prcomp(X)
  expect_equal(
    W,
    sweep(pca$rotation[, 1:L], 2, pca$sdev[1:L], "*")
  )
})

test_that("ARD precision initialises", {
  L <- 4
  alpha <- rgamma(L, shape = 1, rate = 1)

  checkmate::expect_numeric(
    initialise_loading_ard_precision(
      L = L,
      ard_shape = 2, ard_rate = 2,
      ard_precision = NULL,
      perform_checks = TRUE
    ),
    len = L, lower = 0, any.missing = FALSE
  )

  expect_equal(
    alpha,
    initialise_loading_ard_precision(
      L = L,
      ard_shape = 2, ard_rate = 2,
      ard_precision = alpha,
      perform_checks = TRUE
    )
  )
})

test_that("Conditional Correlation", {
  ind_at <- list(
    ord = 1L, nom = 1 + 1:3, con = 5L, fvt = 5 + 1:32
  )
  D_p <- sum(sapply(ind_at, length))

  C_w <- diag(D_p)
  x <- seq(0, 1, length.out = length(ind_at[[4]]))
  d <- abs(outer(x, x, "-"))
  ell <- 1 / (2 * pi)
  C_w[ind_at[[4]], ind_at[[4]]] <-  (exp_quad_kernel(d, 1, ell) + (1e-6 * diag(length(ind_at[[4]])))) / (1 + 1e-6)

  c_star <- compute_scaled_conditional_row_variance_vector(C_w)
  i <- 10
  for(i in 1:D_p) {
    expect_equal(
      c_star[i],
      c(C_w[i, i] - (C_w[i, -i] %*% (solve(C_w[-i, -i] )%*% C_w[-i, i]))),
      tolerance = 1e-4
    )
  }
})

test_that("Conditional mean weights",{
  ind_at <- list(
    ord = 1L, nom = 1 + 1:3, con = 5L, fvt = 5 + 1:32
  )
  D_p <- sum(sapply(ind_at, length))

  C_w <- diag(D_p)
  x <- seq(0, 1, length.out = length(ind_at[[4]]))
  d <- abs(outer(x, x, "-"))
  ell <- 1 / (2 * pi)
  C_w[ind_at[[4]], ind_at[[4]]] <-  (exp_quad_kernel(d, 1, ell) + (1e-6 * diag(length(ind_at[[4]])))) / (1 + 1e-6)

  loading_row_conditional_mean_weight <- compute_loading_row_conditional_mean_weight_matrix(C_w)

  for(i in 1:D_p) {
    expect_equal(
      loading_row_conditional_mean_weight[i, ],
      c(C_w[i, -i] %*% solve(C_w[-i, -i])),
      tolerance = 1e-4
    )
  }
})

test_that("Loading row precision", {
  ind_at <- list(
    ord = 1L, nom = 1 + 1:3, con = 5L, fvt = 5 + 1:32
  )
  D_p <- sum(sapply(ind_at, length))

  C_w <- diag(D_p)
  x <- seq(0, 1, length.out = length(ind_at[[4]]))
  d <- abs(outer(x, x, "-"))
  ell <- 1 / (2 * pi)
  C_w[ind_at[[4]], ind_at[[4]]] <-  (exp_quad_kernel(d, 1, ell) + (1e-6 * diag(length(ind_at[[4]])))) / (1 + 1e-6)

  c_star <- compute_scaled_conditional_row_variance_vector(C_w)

  L <- 4
  alpha <- initialise_loading_ard_precision(
    L = L,
    ard_shape = 2, ard_rate = 2
  )

  P <- 4
  tn <- c("ord", "nom", "con", "fvt")
  tt <- factor(tn, levels = c("ord", "nom", "con", "fvt"))

  lambda <- initialise_precision(
    n_traits = P, trait_names = tn,  trait_type = tt,
    precision_prior_shape = 1, precision_prior_rate = 0.01,
    precision = NULL,
    perform_checks = TRUE
  )
  lambda_vector <- map_precision_to_auxiliary_traits(
    precision = lambda, auxiliary_trait_index = ind_at
  )

  N <- 500
  Z <- matrix(rnorm(N*L), nrow = N, ncol = L)
  ZTZ <- t(Z) %*% Z

  for(i in 1:D_p) {
    expect_equal(compute_loading_row_precision(
      total_individual_specific_latent_trait_outer_product_expectation = ZTZ,
      precision = lambda_vector[i],
      ard_precision = alpha,
      scaled_conditional_row_variance = c_star[i],
      perform_checks = TRUE
    ),
    (lambda_vector[i] * ZTZ)+ (diag(alpha) / c_star[i])
    )
  }
})

test_that("Loading row precision list", {
  ind_at <- list(
    ord = 1L, nom = 1 + 1:3, con = 5L, fvt = 5 + 1:32
  )
  D_p <- sum(sapply(ind_at, length))

  C_w <- diag(D_p)
  x <- seq(0, 1, length.out = length(ind_at[[4]]))
  d <- abs(outer(x, x, "-"))
  ell <- 1 / (2 * pi)
  C_w[ind_at[[4]], ind_at[[4]]] <-  (exp_quad_kernel(d, 1, ell) + (1e-6 * diag(length(ind_at[[4]])))) / (1 + 1e-6)

  c_star <- compute_scaled_conditional_row_variance_vector(C_w)

  L <- 4
  alpha <- initialise_loading_ard_precision(
    L = L,
    ard_shape = 2, ard_rate = 2
  )

  P <- 4
  tn <- c("ord", "nom", "con", "fvt")
  tt <- factor(tn, levels = c("ord", "nom", "con", "fvt"))

  lambda <- initialise_precision(
    n_traits = P, trait_names = tn,  trait_type = tt,
    precision_prior_shape = 1, precision_prior_rate = 0.01,
    precision = NULL,
    perform_checks = TRUE
  )
  lambda_vector <- map_precision_to_auxiliary_traits(
    precision = lambda, auxiliary_trait_index = ind_at
  )

  N <- 500
  Z <- matrix(rnorm(N*L), nrow = N, ncol = L)
  ZTZ <- t(Z) %*% Z

  for(i in 1:D_p) {
    expect_equal(compute_loading_row_precision(
      total_individual_specific_latent_trait_outer_product_expectation = ZTZ,
      precision = lambda_vector[i],
      ard_precision = alpha,
      scaled_conditional_row_variance = c_star[i],
      perform_checks = TRUE
    ),
    (lambda_vector[i] * ZTZ) + (diag(alpha) / c_star[i])
    )
  }

  wwT <- compute_loading_row_precision_list(
    total_individual_specific_latent_trait_outer_product_expectation = ZTZ,
    precision_vector = lambda_vector,
    ard_precision = alpha,
    scaled_conditional_row_variance_vector = c_star,
    perform_checks = TRUE
  )

  for(i in 1:D_p) {
    expect_equal(compute_loading_row_precision(
      total_individual_specific_latent_trait_outer_product_expectation = ZTZ,
      precision = lambda_vector[i],
      ard_precision = alpha,
      scaled_conditional_row_variance = c_star[i],
      perform_checks = TRUE
    ),
    wwT[[i]]
    )
  }

})

test_that("Loading expectation", {
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
    precision = NULL,
    ard_precision = NULL,
    ard_shape = 1, ard_rate = 1,
    loading = NULL, method = "random",
    within_taxon_amplitude = NULL,
    heritable_amplitude = NULL,
    length_scale = 2,
    perform_checks = TRUE
  )

  W_tmp <- plvm$loading_expectation

  W_up <- compute_loading_expectation(
    current_loading_expectation = plvm$loading_expectation,
    loading_row_precision = plvm$loading_row_precision,
    auxiliary_trait = plvm$auxiliary_traits,
    latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
    precision_vector = plvm$precision_vector,
    ard_precision = plvm$ard_precision,
    scaled_conditional_row_variance_vector = plvm$scaled_conditional_loading_row_variance_vector,
    loading_row_conditional_mean_weight = plvm$loading_row_conditional_mean_weight,
    perform_checks = TRUE
  )
  i <- 1
  for (i in 1:D_p) {
    w_star <- c(C_w[i, -i] %*% solve(C_w[-i, -i]) %*% W_tmp[-i, ])
    lxz <- c(plvm$precision_vector[i] *
      t(plvm$auxiliary_traits[, i]) %*%
      plvm$individual_specific_latent_trait_expectation)
    W_tmp[i, ] <- plvm$loading_row_covariance[,, i] %*% (
      lxz +
        (plvm$scaled_conditional_loading_row_variance_vector[i] *
           plvm$ard_precision + w_star)
    )
  }

  expect_equal(
    W_tmp, W_tmp
  )
})

test_that("Loading elbo computation", {
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
    precision = NULL,
    ard_precision = NULL,
    ard_shape = 1, ard_rate = 1,
    loading = NULL, method = "random",
    within_taxon_amplitude = NULL,
    heritable_amplitude = NULL,
    length_scale = 2,
    perform_checks = TRUE
  )

  elbo <- compute_loading_elbo(
    loading_expectation = plvm$loading_expectation,
    loading_row_covariance = plvm$loading_row_covariance,
    ard_precision = plvm$ard_precision,
    loading_prior_correlation_log_det = NULL,
    inv_loading_prior_correlation = NULL,
    loading_prior_correlation = plvm$loading_prior_correlation,
    perform_checks = TRUE
  )

  expect_equal(
    elbo,
    - 0.5 * L * determinant(C_w)$modulus[1] +
      0.5 * D_p * sum(log(plvm$ard_precision)) -
      0.5 * sum(sapply(
        1:L,
        function(i) plvm$ard_precision[i] * sum(diag(
          (plvm$loading_expectation[, i] %*% t(plvm$loading_expectation[, i]) +
             diag(plvm$loading_row_covariance[i, i, ])) %*%
            plvm$inv_loading_prior_correlation
        ))
      )) +
      0.5 * sum(sapply(
        1:D_p,
        function(i){
          determinant(exp(1) * plvm$loading_row_covariance[, , i])$modulus[1]
        }
      ))
  )
})
