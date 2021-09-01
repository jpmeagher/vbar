test_that("initialise auxiliary traits", {
  P <- 4
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
  g <- list(
    ord = function(x){
      ordinal_link(x, cut_off_points = gamma$ord)
    },
    nom = function(x) {
      nominal_link(x, levels = cat$nom)
    },
    con = function(x) x, fvt = function(x) exp(x)
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
        y, mu = matrix(0, N, K[2]),
        n_samples = 1000, return_expectation = FALSE
      )
    },
    con = function(y) data.matrix(y),
    fvt = function(y) log(data.matrix(y))
  )
  gamma <- list(
    ord = c(-Inf, 0, 1, 2, Inf),
    nom = NA, con = NA, fvt = NA
  )
  cat <- list(
    ord = NA,
    nom = factor(levels(mt$nom)),
    con = NA, fvt = NA
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

  for (i in 1:3) {
    expect_equal(
      g[[i]](X[, ind_at[[i]]]), mt[, ind_mt[[i]]]
    )
  }
  i <- 4
  expect_equal(
    g[[i]](X[, ind_at[[i]]]), data.matrix(mt[, ind_mt[[i]]]),
    ignore_attr = TRUE
  )

  expect_equal(
    initialise_auxiliary_traits(
      n_traits = nrow(meta),
      manifest_trait_df = mt,
      trait_names = meta$trait_names,
      trait_type = meta$trait_type,
      trait_levels = meta$trait_levels,
      manifest_trait_index = meta$manifest_trait_index,
      auxiliary_trait_index = meta$auxiliary_trait_index,
      inverse_link_functions = meta$inverse_link_functions,
      auxiliary_traits = X,
      perform_checks = TRUE
    ),
    X
  )
})

test_that("plvm initialised for cavi", {
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

  plvm <- initialise_plvm(
    manifest_trait_df = mt, metadata = meta,
    latent_dimension = L,
    loading_prior_correlation = diag(D_p),
    auxiliary_traits = NULL,
    perform_checks = TRUE
  )

  checkmate::expect_list(plvm)

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

  expect_equal(
    X,
    initialise_plvm(
      manifest_trait_df = mt, metadata = meta,
      latent_dimension = L,
      loading_prior_correlation = diag(D_p),
      auxiliary_traits = X,
      perform_checks = TRUE
    )$auxiliary_traits
  )

})
