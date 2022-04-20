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
#' @param random_seed A single value, interpreted as an integer, or NULL.
#' @param auxiliary_traits An NxD' matrix of real numbers. The initial matrix of
#'   auxiliary traits. Only required when \eqn{method == "pca"}
#' @inheritParams ou_kernel
#'
#' @return A D'xL matrix of real values. The initial loading matrix.
initialise_loading <- function(
  D_prime, L,
  ard_precision, loading_prior_correlation,
  loading = NULL, method = "random", random_seed = NULL,
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
    checkmate::assert_choice(method, choices = c("pca", "random", "varimax"))
  }
  if (!is.null(loading)) return(loading)
  if (method == "random") {
    set.seed(random_seed)
    sigma <- sqrt(1 / ard_precision)
    L_chol <- t(chol(loading_prior_correlation))
    Z <- matrix(stats::rnorm(D_prime * L), nrow = D_prime, ncol = L)
    W_star <- L_chol %*% Z
    W <- sweep(W_star, 2, sigma, "*")
  } else {
    if (is.null(auxiliary_traits)) {
      stop("Auxiliary traits required to initialise PLVM from principal components")
    }
    pca <- stats::prcomp(auxiliary_traits)
    W_star <- pca$rotation[, 1:L]
    W_tmp <- sweep(W_star, 2, pca$sdev[1:L], "*")
    if (method == "pca") {
      W <- W_tmp
    }
    if (method == "varimax") {
      vari <- stats::varimax(W_tmp)
      W <- vari$loadings
    }
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

#' Loading Row Conditional Mean Weight
#'
#' Compute the weight of the remaining loading matrix in the conditional mean
#' for each row of the loading matrix.
#'
#' @inheritParams initialise_loading
#'
#' @return A D'x(D'-1) matrix of positive real numbers. The weight of the
#'   remaining loading matrix in the conditional mean of each row in the loading
#'   matrix.
compute_loading_row_conditional_mean_weight_matrix <- function(
  loading_prior_correlation,
  perform_checks = TRUE
){
  if (perform_checks) {
    checkmate::assert_true(isSymmetric(loading_prior_correlation))
    checkmate::assert_true(matrixcalc::is.positive.definite(loading_prior_correlation))
  }
  C <- loading_prior_correlation
  D_prime <- nrow(C)
  loading_conditional_mean_weights <- sapply(
    1:D_prime,
    function(i){
      C[i, -i] %*% chol2inv(chol(C[-i, -i]))
      }
    )
  t(loading_conditional_mean_weights)
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

#' Compute Loading Expectation
#'
#' Compute the expectation of the PLVM loading under the approximate posterior.
#'
#' @inheritParams compute_loading_row_precision_list
#' @inheritParams compute_auxiliary_trait_elbo
#' @param current_loading_expectation A D'xL matrix of real values. The Loading
#'   Expectation under the previous approximate posterior.
#' @param loading_row_precision A LxLxD' array. The precision of each row of the
#'   loading matrix under the approximate posterior distribution.
#' @param precision_vector A D'-dimensional vector of positive real values. The
#'   precision associated with each auxiliary trait,
#' @param loading_row_conditional_mean_weight A D'x(D'-1) matrix of positive
#'   real numbers. The weight of the remaining loading matrix in the conditional
#'   mean of each row in the loading matrix.
#'
#' @return A D'xL matrix of real values. The Loading expectation under the
#'   approximate posterior distribution.
compute_loading_expectation <- function(
  current_loading_expectation, loading_row_precision,
  auxiliary_traits, latent_trait_expectation,
  precision_vector,
  ard_precision, scaled_conditional_row_variance_vector,
  loading_row_conditional_mean_weight,
  perform_checks = TRUE
){
  N <- nrow(auxiliary_traits)
  D <- nrow(current_loading_expectation)
  L <- ncol(current_loading_expectation)
  if (perform_checks) {
    checkmate::assert_matrix(
      current_loading_expectation, mode = "numeric", any.missing = FALSE
    )
    checkmate::assert_array(
      loading_row_precision, d = 3, mode = "numeric", any.missing = FALSE
    )
    checkmate::assert_set_equal(
      dim(loading_row_precision), c(L, L, D)
    )
    checkmate::assert_matrix(
      auxiliary_traits, ncols = D,
      mode = "numeric", any.missing = FALSE
    )
    checkmate::assert_matrix(
      latent_trait_expectation, nrows = N, ncols = L,
      mode = "numeric", any.missing = FALSE
    )
    checkmate::assert_numeric(
      precision_vector, lower = 0, any.missing = FALSE, len = D
    )
    checkmate::assert_numeric(
      ard_precision, lower = 0, any.missing = FALSE, len = L
    )
    checkmate::assert_matrix(
      loading_row_conditional_mean_weight, nrows = D, ncols = D-1,
      mode = "numeric", any.missing = FALSE
    )
    checkmate::assert_numeric(
      scaled_conditional_row_variance_vector, lower = 0, any.missing = FALSE, len = D
    )
  }
  W <- current_loading_expectation
  LXtZ <- t(auxiliary_traits %*% diag(precision_vector)) %*% latent_trait_expectation
  for (i in 1:D) {
    w_star <- c(t(loading_row_conditional_mean_weight[i, ]) %*% W[-i, ])
    W[i, ] <- solve(
      loading_row_precision[, , i],
      LXtZ[i, ] +  (ard_precision / scaled_conditional_row_variance_vector[i]) * w_star
    )
  }
  W
}

#' Loading ELBO
#'
#' Compute the contribution of the loading matrix Evidence Lower Bound (ELBO) of
#' a PLVM given the approximate posterior distribution for loadings, the prior
#' correlation matrix, and ARD precision hyperparameters.
#'
#' @inheritParams compute_auxiliary_trait_elbo
#' @param loading_row_covariance A LxLxD' array. The covariance of each row of
#'   the loading matrix under the approximate posterior distribution.
#' @param ard_precision An L-dimensional vector of real numbers. The ARD
#'   precision hyper-parameters on the columns of the loading matrix.
#' @param loading_prior_correlation_log_det A real valued scalar. The log
#'   determinant of the loading prior correlation matrix.
#' @param inv_loading_prior_correlation A D'xD' matrix. The inverse of the
#'   loading prior correlation matrix.
#' @param loading_prior_correlation A D'xD' matrix. The loading prior
#'   correlation matrix.
#'
#' @return A real valued scalar. The contribution of the loading matrix to the
#'   ELBO under the approximate posterior distribution.
compute_loading_elbo <- function(
  loading_expectation, loading_row_covariance,
  ard_precision,
  loading_prior_correlation_log_det = NULL,
  inv_loading_prior_correlation = NULL,
  loading_prior_correlation = NULL,
  perform_checks = TRUE
){
  D <- nrow(loading_expectation)
  L <- ncol(loading_expectation)
  if (perform_checks) {
    checkmate::assert_matrix(loading_expectation, mode = "numeric", any.missing = FALSE)
    checkmate::assert_array(loading_row_covariance, mode = "numeric", d = 3, any.missing = FALSE)
    checkmate::assert_set_equal(dim(loading_row_covariance), c(L, L, D))
    checkmate::assert_numeric(ard_precision, len = L, any.missing = FALSE, lower = 0)
    checkmate::assert_matrix(loading_prior_correlation, nrows = D, ncols = D, any.missing = FALSE, null.ok = TRUE)
    checkmate::assert_matrix(inv_loading_prior_correlation, nrows = D, ncols = D, any.missing = FALSE, null.ok = TRUE)
    checkmate::assert_number(loading_prior_correlation_log_det, null.ok = TRUE)
  }
  W_col_outer <- simplify2array(lapply(
    1:L, function(i){
      (loading_expectation[, i] %*% t(loading_expectation[, i])) +
        diag(loading_row_covariance[i, i, ])
      }
  ))
  if (is.null(loading_prior_correlation_log_det)) {
    loading_prior_correlation_log_det <- 2 * sum(log(diag(chol(loading_prior_correlation))))
  }
  if (is.null(inv_loading_prior_correlation)) {
    inv_loading_prior_correlation <- chol2inv(chol(loading_prior_correlation))
  }
  A1 <- - 0.5 * L * loading_prior_correlation_log_det
  A2 <-  0.5 * D * sum(log(ard_precision))
  A3 <- - sum(sapply(
    1:L, function(i){
      0.5 * ard_precision[i] * sum(diag(W_col_outer[, , i] %*% inv_loading_prior_correlation))
    }
  ))
  A4 <- 0.5 * sum(sapply(
    1:D, function(i) 2 * sum(log(diag(chol(exp(1) * loading_row_covariance[, , i]))))
  ))
  A1 + A2 + A3 + A4
}

