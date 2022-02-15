cavi_plvm <- function(
  plvm_list,
  tol = 1e-6, max_iter = 1000,
  n_samples = 1000, random_seed = NULL,
  progress_bar = TRUE,
  perform_checks = TRUE
){
  mod <- plvm_list
  elbo <- c(
    -Inf,
    compute_plvm_elbo(
      mod,
      n_samples = n_samples, random_seed = random_seed,
      perform_checks = perform_checks
    )
  )
  i <- 1
  if (progress_bar) pb <- utils::txtProgressBar(min = 0, max = max_iter, style = 3)
  while (i <= max_iter) {
    #  Phylogenetic GP
    stop("Need a function to update phylogenetic GP hyperparameters")
    # auxiliary traits
    mod$auxiliary_traits <-  update_discrete_auxiliary_traits(
      manifest_trait_df = mod$manifest_trait_df,
      metadata = mod$metadata,
      auxiliary_traits = mod$auxiliary_traits,
      loading_expectation = mod$loading_expectation,
      latent_trait_expectation = mod$individual_specific_latent_trait_expectation,
      n_samples = n_samples, random_seed = random_seed,
      perform_checks = FALSE
    )
    # Precision
    mod$precision <- update_precision(
      precision = mod$precision,
      metadata = mod$metadata,
      auxiliary_traits = mod$auxiliary_traits,
      loading_expectation = mod$loading_expectation,
      latent_trait_expectation = mod$individual_specific_latent_trait_expectation,
      loading_outer_expectation = mod$loading_row_outer_product_expectation,
      latent_trait_outer_expectation = mod$individual_specific_latent_trait_outer_product_expectation,
      perform_checks = FALSE
    )
    mod$precision_vector <- map_precision_to_auxiliary_traits(
      precision = mod$precision,
      auxiliary_trait_index = od$metadata$auxiliary_trait_index,
      perform_checks = FALSE
    )
    # ELBO
    elbo <- c(
      elbo,
      compute_plvm_elbo(
        mod,
        n_samples = n_samples, random_seed = random_seed,
        perform_checks = FALSE
      )
    )
    if (progress_bar) utils::setTxtProgressBar(pb, i)
    i <- i + 1
  }
  if (progress_bar) close(pb)

  list(
    model = mod,
    elbo = elbo
  )
}
