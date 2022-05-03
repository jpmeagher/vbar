#' Ancestral Trait Reconstruction
#'
#' Reconstruct an ancestral trait distribution under the PLVM. Reconstructs the
#' ancestral traits at a specific node on the phylogeny.
#'
#' @param taxon_specific_latent_trait_expectation An L-dimensional vector of
#'   real values. The expected value of taxon-specific latent traits under the
#'   variational distribution.
#' @param taxon_specific_latent_trait_covariance An L-dimensional vector of
#'   positive real values. The variance associated with each taxon-specific
#'   latent trait under the variational distribution.
#' @param loading_expectation A D x L matrix of real values. The expected value
#'   of the Loading matrix under the vaiational distribution.
#' @param precision_vector A D dimensional vector of positive real values. The
#'   independent precision associated with each auxiliary trait.
#' @inheritParams initialise_plvm
#' @param n_samples A positive integer. The number of samples drawn to construct
#'   a Monte-Carlo approximation to the approximate posterior distribution of
#'   nominal traits.
#'
#' @return A named list of length P. For continuous- or function-valued manifest
#'   traits, returns the expectation and marginal standard deviation for each
#'   element of the trait. For ordinal or nomianl traits, returns the
#'   probability associated with each category.
#' @export
variational_ancestral_reconstruction <- function(
  taxon_specific_latent_trait_expectation,
  taxon_specific_latent_trait_covariance,
  loading_expectation,
  precision_vector,
  metadata,
  n_samples = 1000,
  perform_checks = TRUE
  ){
  P <- nrow(metadata)
  L <- length(taxon_specific_latent_trait_expectation)
  D <- length(precision_vector)
  if (perform_checks) {
    checkmate::assert_numeric(taxon_specific_latent_trait_expectation, any.missing = FALSE)
    checkmate::assert_numeric(
      taxon_specific_latent_trait_covariance, len = L, any.missing = FALSE
    )
    checkmate::assert_matrix(
      loading_expectation,
      nrows = D,
      ncols = L,
      any.missing = FALSE
    )
    checkmate::assert_numeric(precision_vector,
                              lower =  0,
                              any.missing = FALSE)
    checkmate::assert_data_frame(metadata)
    checkmate::assert_integerish(n_samples, lower = 1)
  }
  x_ex <-
    c(loading_expectation %*% diag(taxon_specific_latent_trait_expectation))
  x_var <- diag(1 / precision_vector) +
    (
      loading_expectation %*% diag(taxon_specific_latent_trait_covariance) %*% t(loading_expectation)
    )

  ar <- lapply(
    1:P,
    function(i){
      if (metadata$trait_type[i] %in% c("con", "fvt")) {
        tmp <- list()
        tmp$expectation <- x_ex[metadata$auxiliary_trait_index[[i]]]
        tmp$marginal_sd <- sqrt(diag(x_var))[metadata$auxiliary_trait_index[[i]]]
      }
      if (metadata$trait_type[i] == "ord") {
        K <- length(metadata$cut_off_points[[i]]) - 1
        tmp <- list()
        tmp_mu <- x_ex[metadata$auxiliary_trait_index[[i]]]
        tmp_sigma <- sqrt(diag(x_var))[metadata$auxiliary_trait_index[[i]]]
        tmp$prob <-
          stats::pnorm(metadata$cut_off_points[[i]], mean = tmp_mu, sd = tmp_sigma)[1 + 1:K] -
          stats::pnorm(metadata$cut_off_points[[i]], mean = tmp_mu, sd = tmp_sigma)[1:K]
      }
      if (metadata$trait_type[i] == "nom") {
        tmp <- list()
        tmp_mu <- x_ex[metadata$auxiliary_trait_index[[i]]]
        tmp_sigma <-
        x_var[metadata$auxiliary_trait_index[[i]], metadata$auxiliary_trait_index[[i]]]
        sampled_X <- mvnfast::rmvn(n_samples, mu = tmp_mu, sigma = tmp_sigma)
        sampled_y <- apply(sampled_X, 1, which.max)
        tmp$prob <- tabulate(sampled_y) / n_samples
      }
      tmp
    })
  names(ar) <- metadata$trait_names
  ar
}
