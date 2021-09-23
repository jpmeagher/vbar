test_that("elbo is computed", {
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
    perform_checks = TRUE)

  elbo <- compute_plvm_elbo(
    plvm_list = plvm,
    n_samples = 1000, random_seed = NULL,
    perform_checks = TRUE
  )

  expect_equal(
    elbo,
    compute_auxiliary_trait_elbo(
      manifest_trait_df = mt, metadata = meta,
      auxiliary_traits = plvm$auxiliary_traits,
      loading_expectation = plvm$loading_expectation,
      latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
      precision = plvm$precision,
      loading_outer_expectation = plvm$loading_row_outer_product_expectation,
      latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
      n_samples = 1000, random_seed = NULL,
      perform_checks = TRUE
    ) +
      compute_loading_elbo(
        loading_expectation = plvm$loading_expectation,
        loading_row_covariance = plvm$loading_row_covariance,
        ard_precision = plvm$ard_precision,
        loading_prior_correlation_log_det = NULL,
        inv_loading_prior_correlation = NULL,
        loading_prior_correlation = plvm$loading_prior_correlation,
        perform_checks = TRUE
      ) +
      compute_individual_specific_latent_trait_elbo(
        individual_specific_latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
        taxon_id = mt$taxon_id, phy = ph,
        terminal_taxon_specific_latent_trait_expectation = plvm$taxon_specific_latent_trait_expectation[1:S,],
        individual_specific_latent_trait_covariance = plvm$individual_specific_latent_trait_covariance,
        individual_specific_latent_trait_outer_product_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
        terminal_taxon_latent_trait_outer_product_expectation = plvm$taxon_specific_latent_trait_outer_product_expectation[, , 1:S],
        within_taxon_amplitude = plvm$within_taxon_amplitude,
        perform_checks = TRUE
      ) +
      compute_taxon_specific_latent_trait_elbo(
        taxon_specific_latent_trait_expectation = plvm$taxon_specific_latent_trait_expectation,
        taxon_specific_latent_trait_outer_product_expectation = plvm$taxon_specific_latent_trait_outer_product_expectation,
        taxon_specific_latent_trait_covariance = plvm$taxon_specific_latent_trait_covariance,
        phy = ph,
        phylogenetic_gp = plvm$phylogenetic_GP,
        perform_checks = TRUE
      ),
    tolerance = 0.001
  )

  elbo <- compute_plvm_elbo(
    plvm_list = plvm,
    n_samples = 1000, random_seed = 101,
    perform_checks = TRUE
  )

  expect_equal(
    elbo,
    compute_auxiliary_trait_elbo(
      manifest_trait_df = mt, metadata = meta,
      auxiliary_traits = plvm$auxiliary_traits,
      loading_expectation = plvm$loading_expectation,
      latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
      precision = plvm$precision,
      loading_outer_expectation = plvm$loading_row_outer_product_expectation,
      latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
      n_samples = 1000, random_seed = 101,
      perform_checks = TRUE
    ) +
      compute_loading_elbo(
        loading_expectation = plvm$loading_expectation,
        loading_row_covariance = plvm$loading_row_covariance,
        ard_precision = plvm$ard_precision,
        loading_prior_correlation_log_det = NULL,
        inv_loading_prior_correlation = NULL,
        loading_prior_correlation = plvm$loading_prior_correlation,
        perform_checks = TRUE
      ) +
      compute_individual_specific_latent_trait_elbo(
        individual_specific_latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
        taxon_id = mt$taxon_id, phy = ph,
        terminal_taxon_specific_latent_trait_expectation = plvm$taxon_specific_latent_trait_expectation[1:S,],
        individual_specific_latent_trait_covariance = plvm$individual_specific_latent_trait_covariance,
        individual_specific_latent_trait_outer_product_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
        terminal_taxon_latent_trait_outer_product_expectation = plvm$taxon_specific_latent_trait_outer_product_expectation[, , 1:S],
        within_taxon_amplitude = plvm$within_taxon_amplitude,
        perform_checks = TRUE
      ) +
      compute_taxon_specific_latent_trait_elbo(
        taxon_specific_latent_trait_expectation = plvm$taxon_specific_latent_trait_expectation,
        taxon_specific_latent_trait_outer_product_expectation = plvm$taxon_specific_latent_trait_outer_product_expectation,
        taxon_specific_latent_trait_covariance = plvm$taxon_specific_latent_trait_covariance,
        phy = ph,
        phylogenetic_gp = plvm$phylogenetic_GP,
        perform_checks = TRUE
      )
  )
})
