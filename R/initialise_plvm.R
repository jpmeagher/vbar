#' Initialise Auxiliary Traits
#'
#' Initialise auxiliary traits given manifest trait metadata and
#'
#' @inheritParams specify_manifest_trait_metadata
#' @param auxiliary_traits An NxD' matrix of real numbers. The initial matrix of
#'   auxiliary traits. When non null checks that the auxiliary trait matrix is
#'   of the correct type and size.
#'
#' @return An NxD' matrix of real numbers. The initial matrix of auxiliary
#'   traits.
initialise_auxiliary_traits <- function(
  n_traits,
  manifest_trait_df,
  trait_names, trait_type, trait_levels,
  manifest_trait_index,
  auxiliary_trait_index,
  inverse_link_functions,
  auxiliary_traits = NULL,
  perform_checks = TRUE
){
  N <- nrow(manifest_trait_df)
  D_prime <- sum(sapply(auxiliary_trait_index, length))
  if (perform_checks) {
    checkmate::assert_matrix(
      auxiliary_traits, nrows = N, ncols = D_prime,
      mode = "numeric", any.missing = FALSE, null.ok = TRUE
    )
  }
  if (!is.null(auxiliary_traits)) return(auxiliary_traits)
  auxiliary_traits <- matrix(NA, nrow = N, ncol = D_prime)
  for (i in 1:n_traits) {
    auxiliary_traits[, auxiliary_trait_index[[i]]] <-
      inverse_link_functions[[i]](
        manifest_trait_df[, manifest_trait_index[[i]]]
      )
  }
  auxiliary_traits
}

#' Initialise a PLVM
#'
#' Specify a set of initial values for Variational Inference on the Phylogenetic Latent Variable Model (PLVM).
#'
#' @inheritParams initialise_auxiliary_traits
#' @param metadata A data frame. Contains all the metadata required to map a set of manifest traits to the PLVM auxiliary traits..
#' @param latent_dimension An integer. The dimension of latent traits.
#' @param loading_prior_correlation A D'xD' correlation matrix. The correlation within each loading under the Gaussian prior.
#'
#' @seealso specify_manifest_trait_metadata
#'
#' @return A list of intital values for a PLVM..
#' @export
initialise_plvm <- function(
  manifest_trait_df, metadata, latent_dimension,
  loading_prior_correlation,
  auxiliary_traits = NULL,
  perform_checks = TRUE
){
  metadata <- specify_manifest_trait_metadata(
    n_traits = nrow(metadata),
    trait_names = metadata$trait_names,
    trait_type = metadata$trait_type,
    trait_levels = metadata$trait_levels,
    manifest_trait_index = metadata$manifest_trait_index,
    auxiliary_trait_index = metadata$auxiliary_trait_index,
    link_functions = metadata$link_functions,
    inverse_link_functions = metadata$inverse_link_functions,
    cut_off_points = metadata$cut_off_points,
    categories = metadata$categories,
    manifest_trait_df = manifest_trait_df,
    perform_checks = perform_checks
  )
  P <- nrow(metadata)
  N <- nrow(manifest_trait_df)
  L <- latent_dimension
  D <- sum(sapply(metadata$manifest_trait_index, length))
  D_prime <- sum(sapply(metadata$auxiliary_trait_index, length))
  if (perform_checks) {
    checkmate::assert_count(latent_dimension)
    checkmate::assert_matrix(loading_prior_correlation, nrows = D_prime, ncols = D_prime)
    checkmate::assert_true(isSymmetric(loading_prior_correlation))
    checkmate::assert_true(matrixcalc::is.positive.definite(loading_prior_correlation))
  }
  X <- initialise_auxiliary_traits(
    n_traits = nrow(metadata),
    manifest_trait_df = manifest_trait_df,
    trait_names = metadata$trait_names,
    trait_type = metadata$trait_type,
    trait_levels = metadata$trait_levels,
    manifest_trait_index = metadata$manifest_trait_index,
    auxiliary_trait_index = metadata$auxiliary_trait_index,
    inverse_link_functions = metadata$inverse_link_functions,
    auxiliary_traits = auxiliary_traits,
    perform_checks = perform_checks
  )
  out <- list(
    auxiliary_traits = X
  )
  out
}
