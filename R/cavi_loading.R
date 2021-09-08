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
  if (!is.null(loading)) return(loading)
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
#' @param  ard_shape A positive real-valued scalar. The shape of the Gamma
#'   distribution generating initial ARD precision hyper-parameters
#' @param  ard_rate A positive real-valued scalar. The rate of the Gamma
#'   distribution generating initial ARD precision hyper-parameters
#' @param ard_precision An L-dimensional vector of real numbers. The initial ARD
#'   precision hyper-parameters on the loadings. When non-null checks that
#'   \eqn{ard_precision} is of the correct type and size.
#'
#' @return An L-dimensional vector of real numbers. The initial ARD precision
#'   hyper-parameters on the loadings
initialise_loading_ard_precision <- function(
  L,
   ard_shape = 1,  ard_rate = 1,
  ard_precision = NULL,
  perform_checks = TRUE
){
  if (perform_checks) {
    checkmate::assert_count(L, positive = TRUE)
    checkmate::assert_number( ard_shape, lower = 0)
    checkmate::assert_number( ard_rate, lower = 0)
    checkmate::assert_numeric(
      ard_precision, lower = 0, len = L, any.missing = FALSE,
      null.ok = TRUE
    )
  }
  if (is.null(ard_precision)) {
    ard_precision <- stats::rgamma(L, shape =  ard_shape, rate =  ard_rate)
  }
  ard_precision
}

#' Compute the Scaled Conditional Row Variance of the Loading
#'
#' Computes the variance of the loading in each row conditioned on the loading
#' at every other row before allowing for the ARD precision hyper-parameter.
#'
#' @inheritParams initialise_loading
#'
#' @return A D'-dimensional vector of positive real numbers. The scaled
#'   conditional variance of each row in the loading matrix.
compute_scaled_conditional_row_variance_vector <- function(
  loading_prior_correlation,
  perform_checks = TRUE
){
  if (perform_checks) {
    checkmate::assert_true(isSymmetric(loading_prior_correlation))
    checkmate::assert_true(matrixcalc::is.positive.definite(loading_prior_correlation))
  }
  C <- loading_prior_correlation
  D_prime <- nrow(C)
  c_star <- sapply(1:D_prime, function(i){
    c(C[i, i] - C[i, -i] %*% (chol2inv(chol(C[-i, -i])) %*% C[-i, i]))
  })
  c_star
}

#' Compute precision for a row of the Loading Matrix
#'
#' Computes the Gaussian precision associated with a row of the loading matrix
#' under the approximate posterior within Coordinate Ascent Variational
#' Inference (CAVI) for a Phylogenetic Latent Variable Model (PLVM).
#'
#' @inheritParams initialise_loading
#' @inheritParams map_precision_to_auxiliary_traits
#' @param total_individual_specific_latent_trait_outer_product_expectation An
#'   LxL matrix of real values. The sum of expected outer products for each of
#'   the individual specific latent traits under the approximate posterior
#'   within CAVI for the PLVM.
#' @param scaled_conditional_row_variance A positive real-valued scalar. The
#'   variance of the loading in the row conditioned on the loading at every
#'   other row before allowing for the ARD precision hyper-parameter.
#'
#' @return An LxL matrix of real values. The precision associated with a row of
#'   the loading matrix under the approximate distribution within CAVI for a
#'   PLVM.
compute_loading_row_precision <- function(
  total_individual_specific_latent_trait_outer_product_expectation,
  precision,
  ard_precision,
  scaled_conditional_row_variance,
  perform_checks = TRUE
){
  ZTZ <- total_individual_specific_latent_trait_outer_product_expectation
  c_star <- scaled_conditional_row_variance
  L <- length(ard_precision)
  if (perform_checks) {
    checkmate::assert_matrix(
      ZTZ, any.missing = FALSE, nrows = L, ncols = L
    )
    checkmate::assert_number(
      precision, lower = 0, finite = TRUE
    )
    checkmate::assert_numeric(
      ard_precision, any.missing = FALSE, lower = 0, finite = TRUE
    )
    checkmate::assert_number(
      c_star, lower = 0, finite = TRUE
    )
  }
  (precision * ZTZ) + (diag(ard_precision) / c_star)
}


#' Compute precision for each row of the Loading Matrix
#'
#' Computes the Gaussian precision associated with each row of the loading
#' matrix under the approximate posterior within Coordinate Ascent Variational
#' Inference (CAVI) for a Phylogenetic Latent Variable Model (PLVM). Stores the
#' resulting D' precision matrices as a list.
#'
#' @inheritParams compute_loading_row_precision
#' @param precision_vector A D'-dimensional vector of positive real values. The
#'   precision along each dimension of the auxiliary traits.
#' @param scaled_conditional_row_variance_vector A D'-dimensional vector of
#'   positive real values. The variance of the loading in each row conditioned
#'   on the loading at every other row before allowing for the ARD precision
#'   hyper-parameter.
#'
#' @return A list of length D' containing LxL matrices of real values. The
#'   precision associated with each row of the loading matrix under the
#'   approximate distribution within CAVI for a PLVM.
compute_loading_row_precision_list <- function(
  total_individual_specific_latent_trait_outer_product_expectation,
  precision_vector,
  ard_precision,
  scaled_conditional_row_variance_vector,
  perform_checks = TRUE
){
  D_prime <- length(precision_vector)
  ZTZ <- total_individual_specific_latent_trait_outer_product_expectation
  c_star <- scaled_conditional_row_variance_vector
  L <- length(ard_precision)
  if (perform_checks) {
    checkmate::assert_matrix(
      total_individual_specific_latent_trait_outer_product_expectation,
      any.missing = FALSE, nrows = L, ncols = L
    )
    checkmate::assert_numeric(
      precision_vector, lower = 0, finite = TRUE, any.missing = FALSE,
    )
    checkmate::assert_numeric(
      ard_precision, any.missing = FALSE, lower = 0, finite = TRUE
    )
    checkmate::assert_numeric(
      scaled_conditional_row_variance_vector, lower = 0, finite = TRUE,
      len = D_prime, any.missing = FALSE
    )
  }
  lapply(
    1:D_prime,
    function(i){
      compute_loading_row_precision(
        total_individual_specific_latent_trait_outer_product_expectation = ZTZ,
        precision = precision_vector[i],
        ard_precision = ard_precision,
        scaled_conditional_row_variance = c_star[i],
        perform_checks = FALSE
      )
    }
  )
}

