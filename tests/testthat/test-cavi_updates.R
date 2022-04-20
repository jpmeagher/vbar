test_that("precision updated", {
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
    perform_checks = TRUE)

  lambda_updated <- update_precision(
    precision = plvm$precision,
    metadata = meta,
    auxiliary_traits = plvm$auxiliary_traits,
    loading_expectation = plvm$loading_expectation,
    latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
    loading_outer_expectation = plvm$loading_row_outer_product_expectation,
    latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation
  )

  expect_equal(
    lambda_updated["ord"], 1, ignore_attr = T
  )

  expect_equal(
    lambda_updated["nom"], 1, ignore_attr = T
  )

  ind <-  plvm$metadata$auxiliary_trait_index$con
  a_con <- N * length(ind) / 2
  b_con <- 0
  for (i in 1:N) {
    for (j in seq_along(ind)) {
      b_con <- b_con + plvm$auxiliary_traits[i, ind[j]]^2 -
        2 * plvm$auxiliary_traits[i, ind[j]] * (
          t(plvm$loading_expectation[ind[j], ]) %*% plvm$individual_specific_latent_trait_expectation[i, ]
        ) +
        sum(diag(
          plvm$loading_row_outer_product_expectation[, , ind[j]] %*%
            plvm$individual_specific_latent_trait_outer_product_expectation[, , i]
        ))
    }
  }
  b_con <- b_con / 2
  expect_equal(
    lambda_updated["con"], (a_con - 1) / b_con, ignore_attr = T
  )

  ind <-  plvm$metadata$auxiliary_trait_index$fvt
  a_fvt <- N * length(ind) / 2
  b_fvt <- 0
  for (i in 1:N) {
    for (j in seq_along(ind)) {
      b_fvt <- b_fvt + plvm$auxiliary_traits[i, ind[j]]^2 -
        2 * plvm$auxiliary_traits[i, ind[j]] * (
          t(plvm$loading_expectation[ind[j], ]) %*% plvm$individual_specific_latent_trait_expectation[i, ]
        ) +
        sum(diag(
          plvm$loading_row_outer_product_expectation[, , ind[j]] %*%
            plvm$individual_specific_latent_trait_outer_product_expectation[, , i]
        ))
    }
  }
  b_fvt <- b_fvt / 2
  expect_equal(
    lambda_updated["fvt"], (a_fvt - 1) / b_fvt, ignore_attr = T
  )

  i <- 4
  elbo1 <- compute_auxiliary_trait_elbo(
    manifest_trait_df = plvm$manifest_trait_df,
    metadata = plvm$metadata,
    auxiliary_traits = plvm$auxiliary_traits,
    loading_expectation = plvm$loading_expectation,
    latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
    precision = plvm$precision,
    loading_outer_expectation = plvm$loading_row_outer_product_expectation,
    latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
    n_samples = 1000, random_seed = 101,
    perform_checks = TRUE
  )
  elbo2 <- compute_auxiliary_trait_elbo(
    manifest_trait_df = plvm$manifest_trait_df,
    metadata = plvm$metadata,
    auxiliary_traits = plvm$auxiliary_traits,
    loading_expectation = plvm$loading_expectation,
    latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
    precision = lambda_updated,
    loading_outer_expectation = plvm$loading_row_outer_product_expectation,
    latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
    n_samples = 1000, random_seed = 101,
    perform_checks = TRUE
  )
  expect_true(elbo2 > elbo1)

  i <- 3
  elbo1 <- compute_continuous_auxiliary_trait_elbo(
    auxiliary_trait = plvm$auxiliary_traits[, plvm$metadata$auxiliary_trait_index[[i]]],
    loading_expectation = plvm$loading_expectation[plvm$metadata$auxiliary_trait_index[[i]], ],
    latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
    precision = plvm$precision[i],
    loading_outer_expectation = plvm$loading_row_outer_product_expectation[, , plvm$metadata$auxiliary_trait_index[[i]]],
    latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
    perform_checks = TRUE
  )
  elbo2 <- compute_continuous_auxiliary_trait_elbo(
    auxiliary_trait = plvm$auxiliary_traits[, plvm$metadata$auxiliary_trait_index[[i]]],
    loading_expectation = plvm$loading_expectation[plvm$metadata$auxiliary_trait_index[[i]], ],
    latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
    precision = lambda_updated[i],
    loading_outer_expectation = plvm$loading_row_outer_product_expectation[, , plvm$metadata$auxiliary_trait_index[[i]]],
    latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
    perform_checks = TRUE
  )
  expect_true(elbo2 > elbo1)
})

test_that("within-taxon amplitude updated", {
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
    loading = NULL, method = "pca",
    within_taxon_amplitude = NULL,
    heritable_amplitude = NULL,
    length_scale = 2,
    perform_checks = TRUE)

  i <- 1
  expect_equal(
    within_taxon_amplitude_objective(
      par = plvm$within_taxon_amplitude[i], i = i,
      individual_specific_latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
      taxon_id = plvm$manifest_trait_df[, plvm$id_label], phy = plvm$phy,
      terminal_taxon_specific_latent_trait_expectation = plvm$taxon_specific_latent_trait_expectation[1:S,],
      individual_specific_latent_trait_covariance = plvm$individual_specific_latent_trait_covariance,
      individual_specific_latent_trait_outer_product_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
      terminal_taxon_latent_trait_outer_product_expectation = plvm$taxon_specific_latent_trait_outer_product_expectation[, , 1:S],
      within_taxon_amplitude = plvm$within_taxon_amplitude,
      perform_checks = TRUE
    ),
    -compute_individual_specific_latent_trait_elbo(
      individual_specific_latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
      taxon_id = plvm$manifest_trait_df[, plvm$id_label], phy = plvm$phy,
      terminal_taxon_specific_latent_trait_expectation = plvm$taxon_specific_latent_trait_expectation[1:S,],
      individual_specific_latent_trait_covariance = plvm$individual_specific_latent_trait_covariance,
      individual_specific_latent_trait_outer_product_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
      terminal_taxon_latent_trait_outer_product_expectation = plvm$taxon_specific_latent_trait_outer_product_expectation[, , 1:S],
      within_taxon_amplitude = plvm$within_taxon_amplitude,
      perform_checks = TRUE
    )
  )

  current_obj <- within_taxon_amplitude_objective(
    par = plvm$within_taxon_amplitude[i], i = i,
    individual_specific_latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
    taxon_id = plvm$manifest_trait_df[, plvm$id_label], phy = plvm$phy,
    terminal_taxon_specific_latent_trait_expectation = plvm$taxon_specific_latent_trait_expectation[1:S,],
    individual_specific_latent_trait_covariance = plvm$individual_specific_latent_trait_covariance,
    individual_specific_latent_trait_outer_product_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
    terminal_taxon_latent_trait_outer_product_expectation = plvm$taxon_specific_latent_trait_outer_product_expectation[, , 1:S],
    within_taxon_amplitude = plvm$within_taxon_amplitude,
    perform_checks = TRUE
  )
  wt <- plvm$within_taxon_amplitude

  for (i in 1:L) {
    tst <- optim(
      par = wt[i], fn = within_taxon_amplitude_objective,
      i = i,
      individual_specific_latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
      taxon_id = plvm$manifest_trait_df[, plvm$id_label], phy = plvm$phy,
      terminal_taxon_specific_latent_trait_expectation = plvm$taxon_specific_latent_trait_expectation[1:S,],
      individual_specific_latent_trait_covariance = plvm$individual_specific_latent_trait_covariance,
      individual_specific_latent_trait_outer_product_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
      terminal_taxon_latent_trait_outer_product_expectation = plvm$taxon_specific_latent_trait_outer_product_expectation[, , 1:S],
      within_taxon_amplitude = wt,
      perform_checks = FALSE,
      method = "Brent", lower = 0, upper = 1
    )
    expect_true(
      tst$par > 0 & tst$par < 1
    )
    expect_true(
      tst$value < current_obj
    )
    wt[i] <- tst$par
    current_obj <- tst$value
  }
})

test_that("heritable amplitude updated", {
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
    loading = NULL, method = "pca",
    within_taxon_amplitude = NULL,
    heritable_amplitude = NULL,
    length_scale = 2,
    perform_checks = TRUE)

  i <- 1
  expect_equal(
    heritable_amplitude_objective(
      par = plvm$heritable_amplitude[i], i = i,
      heritable_amplitude = plvm$heritable_amplitude, length_scale = plvm$length_scale,
      taxon_specific_latent_trait_expectation = plvm$taxon_specific_latent_trait_expectation,
      taxon_specific_latent_trait_outer_product_expectation = plvm$taxon_specific_latent_trait_outer_product_expectation,
      taxon_specific_latent_trait_covariance = plvm$taxon_specific_latent_trait_covariance,
      phy = plvm$phy,
      phylogenetic_gp = plvm$phylogenetic_GP,
      perform_checks = TRUE
    ),
    -compute_taxon_specific_latent_trait_elbo(
      taxon_specific_latent_trait_expectation = plvm$taxon_specific_latent_trait_expectation,
      taxon_specific_latent_trait_outer_product_expectation = plvm$taxon_specific_latent_trait_outer_product_expectation,
      taxon_specific_latent_trait_covariance = plvm$taxon_specific_latent_trait_covariance,
      phy = plvm$phy,
      phylogenetic_gp = plvm$phylogenetic_GP,
      perform_checks = TRUE
    )
  )

  current_obj <- heritable_amplitude_objective(
    par = plvm$heritable_amplitude[i], i = i,
    heritable_amplitude = plvm$heritable_amplitude, length_scale = plvm$length_scale,
    taxon_specific_latent_trait_expectation = plvm$taxon_specific_latent_trait_expectation,
    taxon_specific_latent_trait_outer_product_expectation = plvm$taxon_specific_latent_trait_outer_product_expectation,
    taxon_specific_latent_trait_covariance = plvm$taxon_specific_latent_trait_covariance,
    phy = plvm$phy,
    phylogenetic_gp = plvm$phylogenetic_GP,
    perform_checks = TRUE
  )
  ha <- plvm$heritable_amplitude
  pgp <- plvm$phylogenetic_GP
  for (i in 1:L) {
    tst <- optim(
      par = ha[i], fn = heritable_amplitude_objective,
      i = i,
      heritable_amplitude = ha, length_scale = plvm$length_scale,
      taxon_specific_latent_trait_expectation = plvm$taxon_specific_latent_trait_expectation,
      taxon_specific_latent_trait_outer_product_expectation = plvm$taxon_specific_latent_trait_outer_product_expectation,
      taxon_specific_latent_trait_covariance = plvm$taxon_specific_latent_trait_covariance,
      phy = plvm$phy,
      phylogenetic_gp = pgp,
      perform_checks = FALSE,
      method = "Brent", lower = 0, upper = 1
    )
    expect_true(
      tst$par > 0 & tst$par < 1
    )
    expect_true(tst$value < current_obj)
    ha[i] <- tst$par
    pgp[, , i] <- reparameterise_phylogenetic_ou(
      phy = plvm$phy,
      heritable_amplitude = ha[i],
      length_scale = plvm$length_scale,
      environmental_amplitude = sqrt(1 - ha[i]^2),
      perform_checks = FALSE
    )
    current_obj <- tst$value
  }
})

test_that("ard loading precision updated", {
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
    perform_checks = TRUE)

  tmp_b <- sapply(
    1:L, function(i){
      1 + (sum(diag(plvm$inv_loading_prior_correlation %*% plvm$loading_col_outer_product_expectation[, , i])) / 2)
    }
  )

  expect_equal(
    (D_p  / 2) / tmp_b,
    update_loading_ard_precision(
      loading_col_outer_product_expectation = plvm$loading_col_outer_product_expectation,
      inv_loading_prior_correlation = plvm$inv_loading_prior_correlation,
      ard_shape = 1,  ard_rate = 1
    )
  )
})

test_that("ordinal trait cut off points updated", {
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
    loading = NULL, method = "pca",
    within_taxon_amplitude = NULL,
    heritable_amplitude = NULL,
    length_scale = 2,
    perform_checks = TRUE)

  i <- 3
  expect_equal(
    ordinal_trait_cut_off_objective(
      par = plvm$metadata$cut_off_points$ord[i], i = i,
      y = plvm$manifest_trait_df$ord,
      cut_off_points = plvm$metadata$cut_off_points$ord,
      loading_expectation = plvm$loading_expectation[plvm$metadata$auxiliary_trait_index$ord, ],
      latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
      loading_outer_expectation = plvm$loading_row_outer_product_expectation[, , plvm$metadata$auxiliary_trait_index$ord],
      latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
      perform_checks = FALSE
    ),
    -compute_ordinal_auxiliary_trait_elbo(
      y = plvm$manifest_trait_df$ord,
      cut_off_points = plvm$metadata$cut_off_points$ord,
      loading_expectation = plvm$loading_expectation[plvm$metadata$auxiliary_trait_index$ord, ],
      latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
      loading_outer_expectation = plvm$loading_row_outer_product_expectation[, , plvm$metadata$auxiliary_trait_index$ord],
      latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation
    )
  )

  for (i in 3:plvm$metadata$trait_levels[1]) {
    tst <- optim(
      par = plvm$metadata$cut_off_points$ord[i], fn = ordinal_trait_cut_off_objective,
      i = i,
      y = plvm$manifest_trait_df$ord,
      cut_off_points = plvm$metadata$cut_off_points$ord,
      loading_expectation = plvm$loading_expectation[plvm$metadata$auxiliary_trait_index$ord, ],
      latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
      loading_outer_expectation = plvm$loading_row_outer_product_expectation[, , plvm$metadata$auxiliary_trait_index$ord],
      latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
      perform_checks = FALSE,
      method = "Brent",
      lower = plvm$metadata$cut_off_points$ord[i-1],
      upper = min(plvm$metadata$cut_off_points$ord[i+1], plvm$metadata$cut_off_points$ord[i] + 10)
    )
    expect_true(
      tst$par > plvm$metadata$cut_off_points$ord[i-1] & tst$par < plvm$metadata$cut_off_points$ord[i+1]
    )
  }
})
