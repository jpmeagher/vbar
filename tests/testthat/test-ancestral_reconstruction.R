test_that("ancestral reconstruction behaves", {
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
    loading = NULL, method = "varimax",
    within_taxon_amplitude = NULL,
    heritable_amplitude = NULL,
    length_scale = 2,
    perform_checks = FALSE)

  iter <- 10
  cavi <- cavi_plvm(
    plvm_list = plvm,
    tol = 1e-1, max_iter = iter,
    n_samples = 1000, random_seed = 202,
    progress_bar = FALSE,
    perform_checks = FALSE
  )

  # ## Match truth to CAVI Loadings
  # R_mat <- IMIFA::Procrustes(
  #   synthetic_trait_model_specification$loading,
  #   cavi$model$loading_expectation,
  # )$R %>% round()
  #
  # matplot(synthetic_trait_model_specification$loading)
  # matplot(cavi$model$loading_expectation %*% R_mat)
  #
  # synthetic_trait_model_specification$individual_specific_latent_traits[, 1]
  # (cavi$model$individual_specific_latent_trait_expectation %*% R_mat)[, 1]
  #
  # plot(
  #   synthetic_trait_model_specification$individual_specific_latent_traits[, 3],
  #   (cavi$model$individual_specific_latent_trait_expectation %*% R_mat)[, 3]
  # )
  #
  # ## Check correlation between individual latent traits and true values
  # i <- 3
  # mod <- lm( synthetic_trait_model_specification$individual_specific_latent_traits[, i] ~
  #       (cavi$model$individual_specific_latent_trait_expectation %*% R_mat)[, i] )
  # summary(mod)
  # cor(synthetic_trait_model_specification$individual_specific_latent_traits[, i],
  #       (cavi$model$individual_specific_latent_trait_expectation %*% R_mat)[, i] )
  #
  # ## Check correlation between taxon latent traits and true values
  # i <- 3
  # plot(synthetic_trait_model_specification$taxon_specific_latent_traits[, i],
  #      (cavi$model$taxon_specific_latent_trait_expectation %*% R_mat)[, i])
  # mod <- lm( synthetic_trait_model_specification$taxon_specific_latent_traits[, i] ~
  #              (cavi$model$taxon_specific_latent_trait_expectation %*% R_mat)[, i] )
  # summary(mod)
  # cor(synthetic_trait_model_specification$taxon_specific_latent_traits[, i],
  #     (cavi$model$taxon_specific_latent_trait_expectation %*% R_mat)[, i])
  #
  # ## Loadings and latent traits look pretty good, so let's look at the ancestral reconstruction

  s <- 1
  ar <- variational_ancestral_reconstruction(
    taxon_specific_latent_trait_expectation = cavi$model$taxon_specific_latent_trait_expectation[s, ],
    taxon_specific_latent_trait_covariance = cavi$model$taxon_specific_latent_trait_covariance[s, ],
    loading_expectation = cavi$model$loading_expectation,
    precision_vector = cavi$model$precision_vector,
    metadata = cavi$model$metadata
    )

  expect_equal(sum(ar$ord$prob), 1)
  expect_equal(sum(ar$nom$prob), 1)

  checkmate::expect_list(
    ar, len = L
  )
})
