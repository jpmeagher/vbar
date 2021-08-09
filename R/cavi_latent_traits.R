#' Compute Individual-specific Latent Trait Precision
#'
#' Computes the precision matrix for the approximate posterior distribution of
#' individual-specific latent traits within the PLVM.
#'
#' @param auxiliary_trait_precision_vector A D-dimensional vector of positive
#'   real values The precision parameters associated with each of the D
#'   auxiliary traits.
#' @param loading_outer_product_expectation An \eqn{LxLxD} dimensional array of
#'   real values. The expected outer product of the loading under the
#'   approximate posterior distribution.
#' @param within_taxon_standard_deviation An L-dimensional vector of real values
#'   on the interval \eqn{[0, 1]}. The within-taxon variation parameters for
#'   each of the L latent traits.
#' @inheritParams ou_kernel
#'
#' @return A symmetric LxL matrix of real values. The precision matrix of the
#'   individual-specific latent traits under the approximate posterior
#'   distribution.
#' @export
compute_individual_specific_latent_trait_precision <- function(
  auxiliary_trait_precision_vector,
  loading_outer_product_expectation,
  within_taxon_standard_deviation,
  perform_checks = TRUE
){
  D <- length(auxiliary_trait_precision_vector)
  L <- length(within_taxon_standard_deviation)
  if (perform_checks) {
    checkmate::assert_numeric(
      auxiliary_trait_precision_vector, any.missing = FALSE
    )
    checkmate::assert_array(
      loading_outer_product_expectation, mode = "numeric", any.missing = FALSE,
      d = 3
    )
    checkmate::assert_true(
      all(dim(loading_outer_product_expectation) == c(L, L, D))
    )
    checkmate::assert_numeric(
      within_taxon_standard_deviation, lower = 0, upper = 1, any.missing = FALSE
    )
  }
  apply(
    sweep(loading_outer_product_expectation, 3, auxiliary_trait_precision_vector, "*"),
    c(1, 2), sum
    ) + diag(1 / within_taxon_standard_deviation^2)
}

#' Compute Individual-specific Latent Trait Expectation
#'
#' Computes the expected value of the individual-specific latent trait  within
#' the PLVM under the approximate posterior distribution.
#'
#' @param auxiliary_trait A D-dimensional vector of real values. The
#'   individual-specific auxiliary trait.
#' @param loading A DxL dimensional matrix of real values. The PLVM loading.
#' @param taxon_specific_latent_trait An L-dimensional vector of real values.
#'   The taxon-specific latent trait values.
#' @inheritParams  compute_individual_specific_latent_trait_precision
#' @param individual_specific_latent_trait_precision A symmetric  LxL matrix of
#'   real values. The precision associated with individual-specific latent
#'   traits.
#'
#' @return An L-dimensional vector of real values. The expected value of the
#'   individual-specific latent trait under the approximate posterior
#'   distributiÒon
#' @export
compute_individual_specific_latent_trait_expectation <- function(
  auxiliary_trait,
  loading,
  taxon_specific_latent_trait,
  auxiliary_trait_precision_vector,
  individual_specific_latent_trait_precision,
  within_taxon_standard_deviation,
  perform_checks = TRUE
){
  D <- length(auxiliary_trait_precision_vector)
  L <- length(within_taxon_standard_deviation)
  if (perform_checks) {
    checkmate::assert_numeric(
      auxiliary_trait, any.missing = FALSE, len = D
    )
    checkmate::assert_matrix(
      loading, mode = "numeric",
      nrows = D, ncols = L, any.missing = FALSE
    )
    checkmate::assert_numeric(
      taxon_specific_latent_trait, any.missing = FALSE, len = L
    )
    checkmate::assert_numeric(
      auxiliary_trait_precision_vector, any.missing = FALSE
    )
    checkmate::assert_matrix(
      individual_specific_latent_trait_precision, mode = "numeric",
      nrows = L, ncols = L, any.missing = FALSE
    )
    checkmate::assert_numeric(
      within_taxon_standard_deviation, lower = 0, upper = 1, any.missing = FALSE
    )
  }
  tmp <- t(loading) %*% diag(auxiliary_trait_precision_vector) %*% auxiliary_trait +
    diag(1 / within_taxon_standard_deviation^2) %*% taxon_specific_latent_trait
  c(solve(
    individual_specific_latent_trait_precision,
    tmp
  ))
}

#' Compute Terminal Taxon-specific Latent Trait Precision
#'
#' Computes the precision of taxon-specific latent traits at terminal nodes of
#' the phylogeny under the approximate posterior distribution for the PLVM.
#'
#'
#' @param N A natural number. The number of individuals witin the taxon.
#' @inheritParams  compute_individual_specific_latent_trait_precision
#' @param conditional_standard_deviation An L-dimensional vector of real values
#'   on the interval \eqn{[0, 1]}. The standard deviation of the taxon-specific
#'   latent trait conditional on its parent taxon.
#'
#' @return An L-dimensional vector of real values. Precision parameters for the
#'   taxo-specific latent traits.
#' @export
compute_terminal_taxon_specific_latent_trait_precision <- function(
  N,
  within_taxon_standard_deviation,
  conditional_standard_deviation,
  perform_checks = TRUE
){
  L <- length(within_taxon_standard_deviation)
  if (perform_checks) {
    checkmate::assert_numeric(
      within_taxon_standard_deviation, lower = 0, upper = 1, any.missing = FALSE
    )
    checkmate::assert_numeric(
      conditional_standard_deviation, lower = 0, upper = 1, any.missing = FALSE
    )
  }
  N / within_taxon_standard_deviation^2 +
    1 / conditional_standard_deviation^2
}
