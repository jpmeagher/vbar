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
})
