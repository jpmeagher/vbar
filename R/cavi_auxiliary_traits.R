#' Ordinal Trait Link function
#'
#' Map an auxiliary trait to ordinal manifest trait.
#'
#' @param x An N-dimensional vector of real numbers. An auxiliary trait.
#' @param cut_off_points A K+1 dimensional ordered vector of values. The cut off
#'   points separating auxiliary traits into ordinal state.
#' @inheritParams ou_kernel
#'
#' @return An N-dimensional vector of ordered factors with K levels. An ordinal
#'   manifest trait.
#' @export
ordinal_link <- function(
  x, cut_off_points,
  perform_checks = TRUE
){
  if (perform_checks) {
    checkmate::assert_numeric(x, any.missing = FALSE)
    checkmate::assert_numeric(cut_off_points, any.missing = FALSE, sorted = TRUE)
  }
  y <- sapply(x, function(z) sum(z > cut_off_points))
  factor(y, levels = 1:(length(cut_off_points) - 1), ordered = TRUE)
}

#' Ordinal Trait Inverse Link Function
#'
#' Map an ordinal manifest trait to an auxiliary trait.
#'
#' @param y An N-dimensional vector of ordered factors with K levels. An ordinal
#'   manifest trait.
#' @inheritParams ordinal_link
#' @param mu An N-dimensional vector of real numbers. The expected value of the
#'   underlying Gaussian random variable with variance 1.
#' @param return_expectation Logical. Should manifest variables map to the
#'   expected auxiliary variable or to values sampled from the truncated
#'   Gaussian distribution.
#' @param random_seed A single value, interpreted as an integer, or NULL.
#'
#' @return An N-dimensional vector of real numbers. An auxiliary trait.
#' @export
ordinal_inverse_link <- function(
  y, cut_off_points,
  mu,
  return_expectation = TRUE, random_seed = NULL,
  perform_checks = TRUE
){
  N <- length(y)
  K <- length(cut_off_points) - 1
  if (perform_checks) {
    checkmate::assert_factor(
      y, n.levels = K, ordered = TRUE, any.missing = FALSE
    )
    checkmate::assert_numeric(cut_off_points, any.missing = FALSE, sorted = TRUE)
    checkmate::assert_numeric(mu, len = N, any.missing = FALSE)
    checkmate::assert_logical(return_expectation)
    checkmate::assert_number(random_seed, null.ok = TRUE)
  }
  y <- as.numeric(y)
  if (return_expectation) {
    x <- RcppTN::etn(
      .mean = mu,
      .low = cut_off_points[y], .high = cut_off_points[y + 1]
    )
  } else {
    set.seed(random_seed)
    x <- RcppTN::rtn(
      .mean = mu,
      .low = cut_off_points[y], .high = cut_off_points[y + 1]
    )
  }
  x
}

#' Nominal Link Function
#'
#' Map the K-dimensional auxiliary trait to a nominal manifest trait.
#'
#' @param X A NxK matrix of real values
#' @param levels An unordered factor vector of length K. The unordered
#'   categories to which manifest traits can belong.
#' @inheritParams ou_kernel
#'
#' @return An N-dimensional vector of unordered factors with K levels. A nominal
#'   manifest trait.
#' @export
nominal_link <- function(
  X, levels,
  perform_checks = TRUE
){
  if (perform_checks) {
    checkmate::assert_matrix(X, mode = "numeric", any.missing = FALSE)
    checkmate::assert_factor(levels, empty.levels.ok = FALSE, n.levels = ncol(X))
  }
  y <- apply(X, 1, which.max)
  levels[y]
}

#' Nominal Probit Normalising Constant
#'
#' Computes a Monte Carlo estimate for probability that a K-dimensional
#' multinomial belongs to category \eqn{i} under the probit model given the mean
#' of the underlying Gaussian distribution with variance 1.
#'
#' @param i A positive integer. The index identifying the multinomial category.
#' @param mu A K-dimensional vector of real numbers. The expected value of the
#'   auxiliary Gaussian mapping to the multinomial.
#' @param n_samples A positive integer. The number of independent samples drawn
#'   to obtain the Monte Carlo estimate.
#' @inheritParams ordinal_inverse_link
#'
#' @return A real number on the unit interval.
#' @export
nominal_probit_normalising_constant <- function(
  i, mu, n_samples = 1000,
  random_seed = NULL,
  perform_checks = TRUE
){
  if (perform_checks == TRUE) {
    K <- length(mu)
    checkmate::assert_integerish(
      i, lower = 1, upper = K, len = 1, any.missing = FALSE
    )
    checkmate::assert_vector(
      mu, strict = TRUE, any.missing = FALSE
    )
    checkmate::assert_count(n_samples, positive = TRUE)
  }
  set.seed(random_seed)
  u <- stats::rnorm(n_samples)
  z <- sapply(u, function(x) x + mu[i] - mu[-i])
  mean(apply(
    stats::pnorm(z), 2, prod
  ))
}

#' Ordinal Probit Normalising Constant
#'
#' Computes the probability that ordinal variables belonging to one of K ordered
#' categories each belong to category \eqn{i} under the probit model given the
#' mean of the underlying Gaussian distribution with variance 1.
#'
#' @inheritParams ordinal_inverse_link
#'
#' @return  A vector of real numbers on the unit interval.
ordinal_probit_normalising_constant <- function(
  y, mu, cut_off_points,
  perform_checks = TRUE
){
  N <- length(mu)
  K <- length(cut_off_points) - 1
  if (perform_checks) {
    checkmate::assert_factor(
      y, n.levels = K, ordered = TRUE, any.missing = FALSE
    )
    checkmate::assert_numeric(mu, len = N, any.missing = FALSE)
    checkmate::assert_numeric(
      cut_off_points, sorted = TRUE, any.missing = FALSE, min.len = 3
    )
    checkmate::assert_set_equal(
      cut_off_points[1:2], c(-Inf, 0)
    )
    checkmate::assert_set_equal(
      cut_off_points[K+1], Inf
    )
  }
  i <- as.numeric(y)
  stats::pnorm(cut_off_points[i+1] - mu) - stats::pnorm(cut_off_points[i] - mu)
}

#' Expectation of multinomial probit auxiliary variables
#'
#' Provides Monte Carlo estimates for the expected value of auxiliary variables in the
#' multinomial probit model given the manifest variable belongs in category \eqn{i} and the K-dimensional mean of the
#' underlying Gaussian distribution with variance 1.
#'
#' @inheritParams nominal_probit_normalising_constant
#'
#' @return A K-dimensional vector of real values.
#' @export
nominal_probit_auxiliary_expectation <- function(
  i, mu, n_samples = 1000,
  random_seed = NULL,
  perform_checks = TRUE
){
  K <- length(mu)
  x_tilde <- rep(NA, K)
  Z <- nominal_probit_normalising_constant(
    i = i, mu = mu, n_samples = n_samples,
    random_seed = random_seed,
    perform_checks = perform_checks
  )
  set.seed(random_seed)
  u <- stats::rnorm(n_samples)
  z <- sapply(u, function(x) x + mu[i] - mu[-i])
  density <- stats::dnorm(z)
  cumulative_density <- stats:: pnorm(z)
  tmp <- sapply(
    1:(K-1),
    function(k){
      mean(apply(rbind(density[k, ], cumulative_density[-k, ]), 2, prod))
    })
  x_tilde[-i] <- mu[-i] - (tmp / Z)
  x_tilde[i] <- mu[i] + sum(mu[-i] - x_tilde[-i])
  x_tilde
}

#' Nominal Inverse Link Function
#'
#' Map a nominal manifest trait to auxiliary traits.
#'
#' @param y An N-dimensional vector of unordered factors with K levels. A nominal
#'   manifest trait.
#' @inheritParams ordinal_inverse_link
#' @inheritParams nominal_probit_auxiliary_expectation
#'
#' @return An NxK matrix of real values.
#' @export
nominal_inverse_link <- function(
  y, mu, n_samples = 1000,
  return_expectation = TRUE, random_seed = NULL,
  perform_checks = TRUE
){
  N <- length(y)
  K <- nlevels(y)
  if (perform_checks) {
    checkmate::assert_factor(y, ordered = FALSE, empty.levels.ok = FALSE, any.missing = FALSE)
    checkmate::assert_matrix(
      mu,
      mode = "numeric", nrows = N, ncols = K,
      any.missing = FALSE
    )
    checkmate::assert_logical(return_expectation)
    checkmate::assert_number(random_seed, null.ok = TRUE)
  }
  y <- as.numeric(y)
  if (return_expectation) {
    X <- t(
      sapply(1:N, function(n){
        nominal_probit_auxiliary_expectation(
          i = y[n], mu = mu[n, ], n_samples = n_samples,
          random_seed = random_seed,
          perform_checks = FALSE
        )
      })
    )
  } else {
    X <- t(
      sapply(1:N, function(i){
        repeat{
          x <- stats::rnorm(K, mean = mu[i, ])
          if (which.max(x) == y[i]) break
        }
        x
      })
    )
  }
  X
}

#' Initialise Auxiliary Traits
#'
#' Initialise auxiliary traits in a PLVM given manifest traits and metadata.
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

#' Initialise Precision
#'
#' Initialise the precision parameters on auxiliary traits within the PLVM.
#'
#' @inheritParams specify_manifest_trait_metadata
#' @param precision_prior_shape A positive real-valued scalar. The shape of the
#'   Gamma prior on the precision.
#' @param precision_prior_rate A positive real-valued scalar. The rate of the
#'   Gamma prior on the precision.
#' @param precision A vector length \eqn{n_traits} taking positive real values.
#'   The Precision with which auxiliary traits are observed. Precision
#'   parameters for discrete traits are fixed at 1. When non-null checks that
#'   the precision parameters have been specified correctly.
#'
#' @return A vector length \eqn{n_traits} taking positive real values. The
#'   Precision with which auxiliary traits are observed. Precision parameters
#'   for discrete traits are fixed at 1.
initialise_precision <- function(
  n_traits, trait_names, trait_type,
  precision_prior_shape = 1, precision_prior_rate = 0.01,
  precision = NULL,
  perform_checks = TRUE
){
  if (perform_checks) {
    checkmate::assert_number(precision_prior_shape, lower = 0)
    checkmate::assert_number(precision_prior_rate, lower = 0)
    checkmate::assert_numeric(
      precision, lower = 0, any.missing = FALSE, null.ok = TRUE,
      names = "named"
      )
    checkmate::assert_numeric(
      precision[trait_type %in% c("ord", "nom")], lower = 1, upper = 1,
      any.missing = FALSE, null.ok = TRUE,
      names = "named"
    )
  }
  if (!is.null(precision)) return(precision)
  precision <- rep(1, n_traits)
  precision[trait_type %in% c("con", "fvt")] <- stats::rgamma(
    sum(trait_type %in% c("con", "fvt")), shape = precision_prior_shape, rate = precision_prior_rate
    )
 names(precision) <- trait_names
 precision
}

#' Map Precision Parameters to Auxiliary Traits
#'
#' Create a vector of precision parameters corresponding to the puxiliary traits.
#'
#' @param precision A vector length \eqn{n_traits} taking positive real values.
#'   The Precision with which auxiliary traits are observed. Precision
#'   parameters for discrete traits are fixed at 1.
#' @inheritParams specify_manifest_trait_metadata
#'
#' @return A D'-dimensional vector of positive real values. The precision along
#'   each dimension of the auxiliary traits.
map_precision_to_auxiliary_traits <- function(
  precision, auxiliary_trait_index,
  perform_checks = TRUE
){
  P <- length(auxiliary_trait_index)
  if (perform_checks) {
    checkmate::assert_numeric(
      precision, lower = 0, any.missing = FALSE, null.ok = TRUE,
      len = P, names = "named"
    )
    checkmate::assert_list(
      auxiliary_trait_index, type = "numeric", any.missing = FALSE,
      names = "named"
    )
    checkmate::assert_set_equal(names(precision), names(auxiliary_trait_index))
  }
  D_prime <- sum(sapply(auxiliary_trait_index, length))
  trait_names <- names(precision)
  precision_vector <- rep(NA, D_prime)
  for (i in 1:P) {
    precision_vector[auxiliary_trait_index[[trait_names[i]]]] <- precision[trait_names[i]]
  }
  precision_vector
}
