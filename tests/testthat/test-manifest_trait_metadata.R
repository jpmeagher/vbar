test_that("manifest trait is recorded as the correct type", {
  P <- 4
  tn <- c("ord", "nom", "con", "fvt")
  tt <- factor(tn, levels = c("ord", "nom", "con", "fvt"))
  ind_mt <- list(
    ord = 2L, nom = 3L, con = 4L, fvt = 5:36
  )
  mt <- cbind(synthetic_traits[, 1:4], fvt = t(simplify2array(synthetic_traits$fvt)))
  K <- c(4, 3, NA, NA)

  i <- 4
  check_trait_type(
    manifest_trait = mt[, ind_mt[[i]]],
    name = tn[i], type = tt[i], levels = K[i]
  )

  expect_error(
    check_trait_type(
      manifest_trait = mt[, ind_mt[[1]]],
      name = tn[1], type = tt[2], trait_levels = K[1]
    )
  )
  expect_error(
    check_trait_type(
      manifest_trait = mt[, ind_mt[[1]]],
      name = tn[1], type = tt[1], trait_levels = K[2]
    )
  )
  expect_error(
    check_trait_type(
      manifest_trait = mt[, ind_mt[[2]]],
      name = tn[2], type = tt[3]
    )
  )
  expect_error(
    check_trait_type(
      manifest_trait = as.numeric(mt[, ind_mt[[2]]]),
      name = tn[2], type = tt[3]
    )
  )
  expect_error(
    check_trait_type(
      manifest_trait = as.numeric(mt[, ind_mt[[2]]]),
      name = tn[2], type = tt[2]
    )
  )
  expect_error(
    check_trait_type(
      manifest_trait = as.character(mt[, ind_mt[[2]]]),
      name = tn[2], type = tt[2]
    )
  )
  expect_error(
    check_trait_type(
      manifest_trait = mt[, ind_mt[[3]]],
      name = tn[3], type = tt[2]
    )
  )
  expect_error(
    check_trait_type(
      manifest_trait = mt[, ind_mt[[4]]],
      name = tn[4], type = tt[3]
    )
  )

})

test_that("manifest trait data frame matches the metadata provided", {
  P <- 4
  tn <- c("ord", "nom", "con", "fvt")
  tt <- factor(tn, levels = c("ord", "nom", "con", "fvt"))
  ind_mt <- list(
    ord = 2L, nom = 3L, con = 4L, fvt = 5:36
  )
  mt <- cbind(synthetic_traits[, 1:4], fvt = t(simplify2array(synthetic_traits$fvt)))
  K <- c(4, 3, NA, NA)

  compare_metadata_to_df(
    n_traits = P, manifest_trait_df = mt,
    trait_names = tn, trait_type = tt,
    trait_levels = K,
    manifest_trait_index = ind_mt
  )

  compare_metadata_to_df(
    n_traits = P, manifest_trait_df = NULL,
    trait_names = tn, trait_type = tt,
    manifest_trait_index = ind_mt
  )

  expect_error(
    compare_metadata_to_df(
      n_traits = P, manifest_trait_df = mt,
      trait_names = c("ord", "nom", "con", "fut"), trait_type = tt,
      trait_levels = K,
      manifest_trait_index = ind_mt
    )
  )
  expect_error(
    compare_metadata_to_df(
      n_traits = P, manifest_trait_df = cbind(mt, fvt.100 = rnorm(nrow(mt))),
      trait_names = c("ord", "nom", "con", "fut"), trait_type = tt,
      trait_levels = K,
      manifest_trait_index = ind_mt
    )
  )
  expect_warning(
    compare_metadata_to_df(
      n_traits = P, manifest_trait_df = cbind(mt, craic = rnorm(nrow(mt))),
      trait_names = tn, trait_type = tt,
      trait_levels = K,
      manifest_trait_index = ind_mt
    )
  )
})

test_that("Manifest trait metadata correctly specified", {
  P <- 4
  tn <- c("ord", "nom", "con", "fvt")
  tt <- factor(tn, levels = c("ord", "nom", "con", "fvt"))
  ind_mt <- list(
    ord = 2L, nom = 3L, con = 4L, fvt = 5:36
  )
  mt <- cbind(synthetic_traits[, 1:4], fvt = t(simplify2array(synthetic_traits$fvt)))
  K <- c(4, 3, NA, NA)
  ind_at <- list(
    ord = 1L, nom = 1 + 1:3, con = 5L, fvt = 5 + 1:32
  )
  g <- list(
    ord = ordinal_link, nom = nominal_link, con = function(x) x, fvt = function(x) exp(x)
  )
  g_inv <- list(
    ord = ordinal_inverse_link, nom = nominal_inverse_link, con = function(y) y, fvt = function(y) log(y)
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


  expect_equal(
    g_inv$fvt(mt[, ind_mt$fvt]),
    log(mt[, ind_mt$fvt])
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
  checkmate::expect_data_frame(meta)

  expect_error(
    specify_manifest_trait_metadata(
      n_traits = P, trait_names = tn, trait_type = tt,
      trait_levels = K,
      manifest_trait_index = ind_mt,
      auxiliary_trait_index = list(
        ord = 1L, nom = 1 + 1:3, con = 5L, fvt = 3 + 1:32
      ),
      link_functions = g,
      inverse_link_functions = g_inv,
      cut_off_points = gamma, categories = cat,
      manifest_trait_df = mt,
      perform_checks = TRUE
    )
  )

  expect_error(
    specify_manifest_trait_metadata(
      n_traits = P, trait_names = tn, trait_type = tt,
      trait_levels = c(4, 3, 1, NA),
      manifest_trait_index = ind_mt, auxiliary_trait_index = ind_at,
      link_functions = g,
      inverse_link_functions = g_inv,
      cut_off_points = gamma, categories = cat,
      manifest_trait_df = mt,
      perform_checks = TRUE
    )
  )

  expect_error(
    specify_manifest_trait_metadata(
      n_traits = P, trait_names = tn, trait_type = tt,
      trait_levels = c(3, 3, NA, NA),
      manifest_trait_index = ind_mt, auxiliary_trait_index = ind_at,
      link_functions = g,
      inverse_link_functions = g_inv,
      cut_off_points = gamma, categories = cat,
      manifest_trait_df = mt,
      perform_checks = TRUE
    )
  )

  expect_error(
    specify_manifest_trait_metadata(
      n_traits = P+1, trait_names = tn, trait_type = tt,
      trait_levels = K,
      manifest_trait_index = ind_mt, auxiliary_trait_index = ind_at,
      link_functions = g,
      inverse_link_functions = g_inv,
      cut_off_points = gamma, categories = cat,
      manifest_trait_df = mt,
      perform_checks = TRUE
    )
  )

  expect_error(
    specify_manifest_trait_metadata(
      n_traits = P, trait_names = tn, trait_type = factor(tn, levels = c("o", "n", "c", "f")),
      trait_levels = K,
      manifest_trait_index = ind_mt, auxiliary_trait_index = ind_at,
      link_functions = g,
      inverse_link_functions = g_inv,
      cut_off_points = gamma, categories = cat,
      manifest_trait_df = mt,
      perform_checks = TRUE
    )
  )

  expect_error(
    specify_manifest_trait_metadata(
      n_traits = P, trait_names = tn, trait_type = factor(tn, levels = c("o", "n", "c", "f")),
      trait_levels = K,
      manifest_trait_index = list(o = 2L, n = 3L, con = 4L, fvt = 5:36),
      auxiliary_trait_index = ind_at,
      link_functions = g,
      inverse_link_functions = g_inv,
      cut_off_points = gamma, categories = cat,
      manifest_trait_df = mt,
      perform_checks = TRUE
    )
  )
})
