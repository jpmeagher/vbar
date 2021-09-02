#' Initialise CAVI for a PLVM
#'
#' Specify a set of initial values for Co-ordinate Ascent Variational Inference
#' (CAVI) on the Phylogenetic Latent Variable Model (PLVM).
#'
#' @inheritParams initialise_auxiliary_traits
#' @inheritParams initialise_loading
#' @inheritParams initialise_loading_ard_precision
#' @param metadata A data frame. Contains all the metadata required to map a set
#'   of manifest traits to the PLVM auxiliary traits..
#'
#' @seealso specify_manifest_trait_metadata
#'
#' @return A list of intital values for a PLVM..
#' @export
initialise_plvm <- function(
  manifest_trait_df, metadata, L,
  loading_prior_correlation,
  auxiliary_traits = NULL,
  ard_precision = NULL,
  ard_prior_shape = 1, ard_prior_rate = 1,
  loading = NULL, method = "random",
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
  D <- sum(sapply(metadata$manifest_trait_index, length))
  D_prime <- sum(sapply(metadata$auxiliary_trait_index, length))
  # Auxiliary Traits
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
  # Loading
  alpha <- initialise_loading_ard_precision(
    L = L,
    ard_prior_shape = ard_prior_shape, ard_prior_rate = ard_prior_rate,
    ard_precision = ard_precision,
    perform_checks = perform_checks
  )
  W <- initialise_loading(
    D_prime = D_prime, L = L,
    ard_precision = alpha,
    loading_prior_correlation = loading_prior_correlation,
    loading = loading, method = method,
    auxiliary_traits = X,
    perform_checks = perform_checks
  )

  out <- list(
    auxiliary_traits = X,
    ard_prior_shape = ard_prior_shape, ard_prior_rate = ard_prior_rate,
    ard_precision = alpha,
    loading_prior_correlation = loading_prior_correlation,
    loading = W
  )
  out
}
