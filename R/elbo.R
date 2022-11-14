#' PLVM ELBO
#'
#' Compute the Evidence Lower Bound (ELBO) for a specific a PLVM.Ê
#'
#' @param plvm_list A list defining the Phylogenetic Latent variable Model
#' @inheritParams compute_auxiliary_trait_elbo
#'
#' @return A scalar. The Evidence Lower Bound for a given Phylogenetic Latent
#'   Variable Model
#' @export
#'
#' @seealso initialise_plvm
compute_plvm_elbo <- function(
  plvm_list,
  n_samples = 1000, random_seed = NULL,
  perform_checks = TRUE
){
  S <- length(plvm_list$phy$tip.label)
  elbo <- compute_auxiliary_trait_elbo(
    manifest_trait_df = plvm_list$manifest_trait_df,
    metadata = plvm_list$metadata,
    auxiliary_traits = plvm_list$auxiliary_traits,
    loading_expectation = plvm_list$loading_expectation,
    latent_trait_expectation = plvm_list$individual_specific_latent_trait_expectation,
    precision = plvm_list$precision,
    loading_outer_expectation = plvm_list$loading_row_outer_product_expectation,
    latent_trait_outer_expectation = plvm_list$individual_specific_latent_trait_outer_product_expectation,
    n_samples = n_samples, random_seed = random_seed,
    perform_checks = perform_checks
  )
  elbo <- elbo +
    compute_loading_elbo(
      loading_expectation = plvm_list$loading_expectation,
      loading_row_covariance = plvm_list$loading_row_covariance,
      ard_precision = plvm_list$ard_precision,
      loading_prior_correlation_log_det = plvm_list$loading_prior_correlation_log_det,
      inv_loading_prior_correlation = plvm_list$inv_loading_prior_correlation,
      loading_prior_correlation = plvm_list$loading_prior_correlation,
      perform_checks = perform_checks
      )
  elbo <- elbo +
    compute_individual_specific_latent_trait_elbo(
      individual_specific_latent_trait_expectation = plvm_list$individual_specific_latent_trait_expectation,
      taxon_id = plvm_list$manifest_trait_df[, plvm_list$id_label], phy = plvm_list$phy,
      terminal_taxon_specific_latent_trait_expectation = plvm_list$taxon_specific_latent_trait_expectation[1:S, , drop = F],
      individual_specific_latent_trait_covariance = plvm_list$individual_specific_latent_trait_covariance,
      individual_specific_latent_trait_outer_product_expectation = plvm_list$individual_specific_latent_trait_outer_product_expectation,
      terminal_taxon_latent_trait_outer_product_expectation = plvm_list$taxon_specific_latent_trait_outer_product_expectation[, , 1:S, drop = F],
      within_taxon_amplitude = plvm_list$within_taxon_amplitude,
      perform_checks = perform_checks
    )
  elbo <- elbo +
    compute_taxon_specific_latent_trait_elbo(
      taxon_specific_latent_trait_expectation = plvm_list$taxon_specific_latent_trait_expectation,
      taxon_specific_latent_trait_outer_product_expectation = plvm_list$taxon_specific_latent_trait_outer_product_expectation,
      taxon_specific_latent_trait_covariance = plvm_list$taxon_specific_latent_trait_covariance,
      phy = plvm_list$phy,
      phylogenetic_gp = plvm_list$phylogenetic_GP,
      perform_checks = perform_checks
    )
  elbo
}
