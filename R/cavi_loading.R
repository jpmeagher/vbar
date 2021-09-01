#' Initialis Loading
#'
#' Initialise loadings for a PLVM. Allows for either a random initialisation,
#' which does not depend on the auxiliary traits, or initialises with the first
#' \eqn{L} principle components, which does.
#'
#' @param D_prime An integer. The dimension of auxiliary traits traits.
#' @param L An integer. The dimension of latent traits.
#' @param loading_prior_correlation A D'xD' correlation matrix. The correlation
#'   within loadings under the Gaussian prior.
#' @param ard_precision An L-dimensional vector of real numbers. The initial ARD
#'   precision hyper-parameters on the loadings.
#' @param loading A D'xL matrix of real values. The initial loading matrix. When
#'   non null checks that the loading matrix is of the correct type and size.
#' @param method A string. Either \eqn{"random"} or \eqn{"pca"} indicating the
#'   type of initialisation adopted.
#' @param auxiliary_traits An NxD' matrix of real numbers. The initial matrix of
#'   auxiliary traits. Only required when \eqn{method == "pca"}
#' @inheritParams ou_kernel
#'
#' @return A D'xL matrix of real values. The initial loading matrix.
initialise_loading <- function(
  D_prime, L,
  ard_precision, loading_prior_correlation,
  loading = NULL, method = "random",
  auxiliary_traits = NULL,
  perform_checks = TRUE
){
  if (perform_checks) {
    checkmate::assert_count(D_prime, positive = TRUE)
    checkmate::assert_count(L)
    checkmate::assert_numeric(
      ard_precision, lower = 0, len = L, any.missing = FALSE
    )
    checkmate::assert_matrix(loading_prior_correlation, nrows = D_prime, ncols = D_prime)
    checkmate::assert_true(isSymmetric(loading_prior_correlation))
    checkmate::assert_true(matrixcalc::is.positive.definite(loading_prior_correlation))
    checkmate::assert_matrix(
      auxiliary_traits, ncols = D_prime,
      mode = "numeric", any.missing = FALSE, null.ok = TRUE
    )
    checkmate::assert_choice(method, choices = c("pca", "random"))
  }
  if (method == "pca") {
    if (is.null(auxiliary_traits)) {
      stop("Auxiliary traits required to initialise PLVM at principal components")
    }
    pca <- stats::prcomp(auxiliary_traits)
    W_star <- pca$rotation[, 1:L]
    W <- sweep(W_star, 2, pca$sdev[1:L], "*")
  }
  if (method == "random") {
    sigma <- sqrt(1 / ard_precision)
    L_chol <- t(chol(loading_prior_correlation))
    Z <- matrix(stats::rnorm(D_prime * L), nrow = D_prime, ncol = L)
    W_star <- L_chol %*% Z
    W <- sweep(W_star, 2, sigma, "*")
  }
  W
}

#' Initialise ARD Precision
#'
#' Initialiase the Automatic Relevance Determination hyper-parameter in the
#' Gaussian process prior on PLVM loadings.
#'
#' @inheritParams initialise_loading
#' @param ard_prior_shape A positive real-valued scalar. The shape of the Gamma
#'   prior on the ARD precision hyper-parameters
#' @param ard_prior_rate A positive real-valued scalar. The rate of the Gamma
#'   prior on the ARD precision hyper-parameters
#' @param ard_precision An L-dimensional vector of real numbers. The initial ARD
#'   precision hyper-parameters on the loadings. When non-null checks that
#'   \eqn{ard_precision} is of the correct type and size.
#'
#' @return An L-dimensional vector of real numbers. The initial ARD precision
#'   hyper-parameters on the loadings
initialise_loading_ard_precision <- function(
  L,
  ard_prior_shape = 1, ard_prior_rate = 1,
  ard_precision = NULL,
  perform_checks = TRUE
){
  if (perform_checks) {
    checkmate::assert_count(L, positive = TRUE)
    checkmate::assert_number(ard_prior_shape, lower = 0)
    checkmate::assert_number(ard_prior_rate, lower = 0)
    checkmate::assert_numeric(
      ard_precision, lower = 0, len = L, any.missing = FALSE,
      null.ok = TRUE
    )
  }
  if (is.null(ard_precision)) {
    ard_precision <- stats::rgamma(L, shape = ard_prior_shape, rate = ard_prior_rate)
  }
  ard_precision
}
