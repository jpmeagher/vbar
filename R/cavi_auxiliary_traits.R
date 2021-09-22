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
#' @param log_out Logical. Should the log probability be returned.
#'
#' @return A (log) probability scalar.
#' @export
nominal_probit_normalising_constant <- function(
  i, mu, n_samples = 1000,
  random_seed = NULL,
  log_out = FALSE,
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
  Z <- mean(apply(
    stats::pnorm(z), 2, prod
  ))
  if (log_out) return(log(Z))
  Z
}

#' Ordinal Probit Normalising Constant
#'
#' Computes the probability that ordinal variables belonging to one of K ordered
#' categories each belong to category \eqn{i} under the probit model given the
#' mean of the underlying Gaussian distribution with variance 1.
#'
#' @inheritParams ordinal_inverse_link
#' @param log_out Logical. Should the log probability be returned.
#'
#' @return  A vector of (log) probabilities.
ordinal_probit_normalising_constant <- function(
  y, mu, cut_off_points,
  log_out = TRUE,
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
  log_Z <- stats::pnorm(cut_off_points[i+1], mean = mu, log.p = TRUE) +
    log1p(-exp(
      stats::pnorm(cut_off_points[i], mean = mu, log.p = TRUE) -
        stats::pnorm(cut_off_points[i+1], mean = mu, log.p = TRUE)
    ))
  if (log_out) return(log_Z)
  exp(log_Z)
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

#' Ordinal Trait ELBO
#'
#' Compute the contribution of an ordinal trait to the Evidence Lower Bound
#' (ELBO) of a PLVM given the approximate posterior distribution for loadings
#' and latent traits.
#'
#' @inheritParams ordinal_probit_normalising_constant
#' @param loading_expectation A L-dimensional vector of real numbers, The row of
#'   the expected loading matrix corresponding to the ordinal trait.
#' @param latent_trait_expectation An NxL matrix of real values. The expected
#'   individual specific latent traits.
#' @param loading_outer_expectation A LxL matrix. The expected outer product for
#'   the row of the expected loading matrix corresponding to the ordinal trait.
#' @param latent_trait_outer_expectation A LxLxN array, The expected outer
#'   product of individual specific latent traits.
#'
#' @return A real valued scalar. The contribution of the ordinal trait to the
#'   ELBO.
compute_ordinal_auxiliary_trait_elbo <- function(
  y, cut_off_points,
  loading_expectation, latent_trait_expectation,
  loading_outer_expectation, latent_trait_outer_expectation,
  perform_checks = TRUE
){
  K <- length(cut_off_points) - 1
  N <- length(y)
  L <- length(loading_expectation)
  if (perform_checks) {
    checkmate::assert_vector(
      loading_expectation, strict = T, any.missing = FALSE
    )
    checkmate::assert_matrix(
      latent_trait_expectation, mode = "numeric",
      nrows = N, ncols = L, any.missing = FALSE
    )
    checkmate::assert_matrix(
      loading_outer_expectation, mode = "numeric",
      any.missing = FALSE, ncols = L, nrows = L
    )
    checkmate::assert_array(
      latent_trait_outer_expectation,
      mode = "numeric", any.missing = FALSE, d = 3
    )
  }
  log_Z <- ordinal_probit_normalising_constant(
    y = y, mu = latent_trait_expectation %*% loading_expectation,
    cut_off_points = cut_off_points, log_out = TRUE,
    perform_checks = perform_checks
  )
  m_ex_2 <- (latent_trait_expectation %*% loading_expectation)^2
  sum_m_2_ex <- sum(diag(loading_outer_expectation %*% apply(latent_trait_outer_expectation, c(1, 2), sum)))
  sum(log_Z) + (0.5 * (sum(m_ex_2) - sum_m_2_ex))
}

#' Nominal Trait ELBO
#'
#' Compute the contribution of a nominal trait to the Evidence Lower Bound
#' (ELBO) of a PLVM given the approximate posterior distribution for loadings
#' and latent traits.
#'
#' @inheritParams nominal_inverse_link
#' @inheritParams nominal_probit_normalising_constant
#' @param loading_expectation A KxL matrix of real numbers, The row of
#'   the expected loading matrix corresponding to the nominal trait.
#' @param latent_trait_expectation An NxL matrix of real values. The expected
#'   individual specific latent traits.
#' @param loading_outer_expectation A LxLxK array. The expected outer products for
#'   the row of the expected loading matrix corresponding to the nominal trait.
#' @param latent_trait_outer_expectation A LxLxN array. The expected outer
#'   product of individual specific latent traits.
#'
#' @return A real valued scalar. The contribution of the nominal trait to the
#'   ELBO.
compute_nominal_auxiliary_trait_elbo <- function(
  y, n_samples = 1000, random_seed = NULL,
  loading_expectation, latent_trait_expectation,
  loading_outer_expectation, latent_trait_outer_expectation,
  perform_checks = TRUE
){
  K <- nrow(loading_expectation)
  L <- ncol(loading_expectation)
  N <- length(y)
  if (perform_checks) {
    checkmate::assert_factor(y, n.levels = K, ordered = FALSE, any.missing = FALSE)
    checkmate::assert_matrix(
      loading_expectation, mode = "numeric", any.missing = FALSE
    )
    checkmate::assert_matrix(
      latent_trait_expectation, mode = "numeric",
      ncols = L, any.missing = FALSE
    )
    checkmate::assert_array(
      loading_outer_expectation, mode = "numeric", d = 3,
      any.missing = FALSE
    )
    checkmate::assert_array(
      latent_trait_outer_expectation,
      mode = "numeric", any.missing = FALSE, d = 3
    )
  }
  num_y <- as.numeric(y)
  m_ex <- latent_trait_expectation %*% t(loading_expectation)
  log_Z <- sapply(
    1:N, function(j){
      if (is.null(random_seed)) {
        rs <- NULL
      } else {
        rs <- random_seed + j
      }
      nominal_probit_normalising_constant(
        num_y[j], mu = m_ex[j, ], n_samples = n_samples,
        random_seed = rs,
        log_out = TRUE,
        perform_checks = (j == 1)
      )
    }
  )
  sum_m_2_ex <- sum(diag(
    apply(loading_outer_expectation, c(1, 2), sum) %*% apply(latent_trait_outer_expectation, c(1, 2), sum)
  ))
  sum(log_Z) + (0.5 * (sum(m_ex^2) - sum_m_2_ex))
}

#' Continuous Trait ELBO
#'
#' Compute the contribution of a continuous trait to the Evidence Lower Bound
#' (ELBO) of a PLVM given the approximate posterior distribution for loadings
#' and latent traits
#'
#' @param auxiliary_trait Either an N-dimensional vector or a NxD' matrix of
#'   real values. The auxiliary trait measurements associated with a single
#'   scalar- or function-valued trait.
#' @param loading_expectation Either a L-dimensional vector or D'xL matrix of
#'   real numbers, The row(s) of the expected loading matrix corresponding to
#'   the continuous trait.
#' @param latent_trait_expectation An NxL matrix of real values. The expected
#'   individual specific latent traits.
#' @param loading_outer_expectation A LxL matrix or LxLxD' array. The expected
#'   outer product(s) for the row(s) of the expected loading matrix
#'   corresponding to the continuous trait.
#' @param latent_trait_outer_expectation A LxLxN array. The expected outer
#'   product of individual specific latent traits.
#' @param precision A positive real valued scalar. The precision with which the
#'   auxiliary trait is observed.
#' @inheritParams ou_kernel
#'
#' @return A scalar. The contribution to the ELBO made by the continuous trait.
compute_continuous_auxiliary_trait_elbo <- function(
  auxiliary_trait,
  loading_expectation, latent_trait_expectation, precision,
  loading_outer_expectation, latent_trait_outer_expectation,
  perform_checks = TRUE
){
  X <- as.matrix(auxiliary_trait)
  N <- nrow(X)
  D <- ncol(X)
  L <- ncol(latent_trait_expectation)
  if (perform_checks) {
    checkmate::assert_matrix(
      X, mode = "numeric", any.missing = FALSE
    )
    checkmate::assert(
      checkmate::check_numeric(
        loading_expectation, any.missing = FALSE, len = L
      ),
      checkmate::check_matrix(
        loading_expectation, mode = "numeric", any.missing = FALSE,
        nrows = D, ncols = L
      )
    )
    checkmate::assert_number(precision, lower = 0)
    checkmate::assert_matrix(
      latent_trait_expectation, mode = "numeric",
      ncols = L, any.missing = FALSE
    )
    checkmate::assert(
      checkmate::check_matrix(
        loading_outer_expectation, mode = "numeric", any.missing = FALSE,
        nrows = L, ncols = L
      ),
      checkmate::check_array(
        loading_outer_expectation, mode = "numeric", d = 3,
        any.missing = FALSE
      )
    )
  }
  if (D == 1) loading_expectation = t(loading_expectation)
  A1 <- - N * D * log(2 * pi) / 2
  A2 <- - N * D * log(precision) / 2
  A3.1 <- sum(diag(t(X) %*% X))
  A3.2 <- - 2 * sum(diag(t(X) %*% (latent_trait_expectation %*% t(loading_expectation))))
  A3.3 <- sum(diag(
    apply(loading_outer_expectation, c(1, 2), sum) %*% apply(latent_trait_outer_expectation, c(1, 2), sum)
  ))
  A3 <- - precision * (A3.1 + A3.2 + A3.3) / 2
  unname(A1 + A2 + A3)
}

#' Auxiliary Trait ELBO
#'
#' Compute the contribution of auxiliary traits to the Evidence Lower Bound
#' (ELBO) of a PLVM given the approximate posterior distribution for loadings
#' and latent traits
#'
#' @inheritParams initialise_plvm
#' @param auxiliary_traits An NxD' matrix of real numbers. The auxiliary traits.
#' @param loading_expectation A D'xL matrix of real numbers, The expected
#'   loading matrix.
#' @param latent_trait_expectation An NxL matrix of real values. The expected
#'   individual specific latent traits.
#' @param loading_outer_expectation A LxLxD' array. The expected outer products
#'   of the expected loading matrix.
#' @param latent_trait_outer_expectation A LxLxN array. The expected outer
#'   product of individual specific latent traits.
#' @param precision A P-dimensional vector of positive real values. The
#'   precision with which auxiliary traits are observed.
#' @inheritParams compute_nominal_auxiliary_trait_elbo
#' @inheritParams compute_ordinal_auxiliary_trait_elbo
#'
#' @return A scalar. The contribution of Auxiliary traits to the ELBO.
compute_auxiliary_trait_elbo <- function(
  manifest_trait_df, metadata,
  auxiliary_traits,
  loading_expectation, latent_trait_expectation,
  precision,
  loading_outer_expectation, latent_trait_outer_expectation,
  n_samples = 1000, random_seed = NULL,
  perform_checks = TRUE
){
  P <- nrow(metadata)
  elbo <- rep(NA, P)
  for (i in 1:P) {
    if (metadata$trait_type[i] %in% c("con", "fvt")) {
      elbo[i] <- compute_continuous_auxiliary_trait_elbo(
        auxiliary_trait = auxiliary_traits[, metadata$auxiliary_trait_index[[i]]],
        loading_expectation = loading_expectation[metadata$auxiliary_trait_index[[i]], ],
        latent_trait_expectation = latent_trait_expectation,
        precision = precision[i],
        loading_outer_expectation = loading_outer_expectation[, , metadata$auxiliary_trait_index[[i]]],
        latent_trait_outer_expectation = latent_trait_outer_expectation,
        perform_checks = perform_checks
      )
    } else {
      if (metadata$trait_type[i] == "nom") {
        elbo[i] <- compute_nominal_auxiliary_trait_elbo(
          y = manifest_trait_df[, metadata$manifest_trait_index[[i]]],
          n_samples = n_samples, random_seed = random_seed,
          loading_expectation = loading_expectation[metadata$auxiliary_trait_index[[i]], ],
          latent_trait_expectation = latent_trait_expectation,
          loading_outer_expectation = loading_outer_expectation[, , metadata$auxiliary_trait_index[[i]]],
          latent_trait_outer_expectation = latent_trait_outer_expectation,
          perform_checks = perform_checks
        )
      }
      if (metadata$trait_type[i] == "ord") {
        elbo[i] <- compute_ordinal_auxiliary_trait_elbo(
          y = manifest_trait_df[, metadata$manifest_trait_index[[i]]],
          cut_off_points = metadata$cut_off_points[[i]],
          loading_expectation = loading_expectation[metadata$auxiliary_trait_index[[i]], ],
          latent_trait_expectation = latent_trait_expectation,
          loading_outer_expectation = loading_outer_expectation[, , metadata$auxiliary_trait_index[[i]]],
          latent_trait_outer_expectation = latent_trait_outer_expectation,
          perform_checks = perform_checks
        )
      }
    }
  }
  sum(elbo)
}
