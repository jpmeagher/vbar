#' Update auxiliary trait precision
#'
#' Update the precision associated with each auxiliary trait within the PLVM by
#' optimising the ELBO.
#'
#' @inheritParams compute_auxiliary_trait_elbo
#'
#' @return A P-dimensional vector. The precision associated with each auxiliary
#'   trait.
update_precision <- function(
  precision,
  metadata,
  auxiliary_traits,
  loading_expectation, latent_trait_expectation,
  loading_outer_expectation, latent_trait_outer_expectation
){
  P <- nrow(metadata)
  N <- nrow(auxiliary_traits)
  X <- auxiliary_traits
  M <- latent_trait_expectation %*% t(loading_expectation)
  lambda <- precision
  for (i in 1:P) {
    if (metadata$trait_type[i] %in% c("con", "fvt")) {
      ind <- metadata$auxiliary_trait_index[[i]]
      D <- length(ind)
      alpha <- N * D / 2
      beta <- (
        sum(diag(t(X[, ind]) %*% X[, ind])) -
          (2 * sum(diag(t(X[, ind]) %*% M[, ind]))) +
          sum(diag(
            apply(loading_outer_expectation[, , ind], c(1, 2), sum) %*%
              apply(latent_trait_outer_expectation, c(1, 2), sum)
          ))
      ) / 2
      lambda[i] <- (alpha - 1) / beta
    }
  }
  lambda
}

#' Objective function for within-taxon amplitude
#'
#' Update the within-taxon amplitude by optimising the ELBO for
#' individual-specific latent traits.
#'
#' @param  par A real value on the unit interval. The within-taxon amplitude.
#' @param i An integer value between 1 and L. The index for the within-taxon
#'   amplitude to be optimised.
#' @inheritParams compute_individual_specific_latent_trait_elbo
#'
#' @return A scalar. The negative ELBO associated with individual specific
#'   latent traits scaled to account for the logit transform of the within-taxon
#'   amplitude.
within_taxon_amplitude_objective <- function(
  par, i,
  individual_specific_latent_trait_expectation,
  taxon_id, phy,
  terminal_taxon_specific_latent_trait_expectation,
  individual_specific_latent_trait_covariance,
  individual_specific_latent_trait_outer_product_expectation,
  terminal_taxon_latent_trait_outer_product_expectation,
  within_taxon_amplitude,
  perform_checks = TRUE
  ){
  wta <- within_taxon_amplitude
  # tmp_par <- stats::plogis(par)
  # wta[i] <- tmp_par
  wta[i] <- par
  obj <- compute_individual_specific_latent_trait_elbo(
    individual_specific_latent_trait_expectation = individual_specific_latent_trait_expectation,
    taxon_id = taxon_id, phy = phy,
    terminal_taxon_specific_latent_trait_expectation = terminal_taxon_specific_latent_trait_expectation,
    individual_specific_latent_trait_covariance = individual_specific_latent_trait_covariance,
    individual_specific_latent_trait_outer_product_expectation = individual_specific_latent_trait_outer_product_expectation,
    terminal_taxon_latent_trait_outer_product_expectation = terminal_taxon_latent_trait_outer_product_expectation,
    within_taxon_amplitude = wta,
    perform_checks = perform_checks
  )
  # - (obj + log(tmp_par) + log(1 - tmp_par))
  - obj
}

#' Objective function for heritable amplitude
#'
#' Update the heritable amplitude by optimising the ELBO for taxon-specific
#' latent traits.
#'
#' @param  par A real value on the unit interval. The heritable amplitude.
#' @param i An integer value between 1 and L. The index for the heritable
#'   amplitude to be optimised.
#' @param heritable_amplitude An L-dimensional vector of values on the unit
#'   interval. The current amplitude of the phylogenetic OU process. Note that
#'   this is the square root of the variance and so should be thought of as a
#'   standard deviation parameter.
#' @param length_scale A positive real-valued scalar. The amplitude of the
#'   phylogenetic OU process. Note that this is the square root of the variance
#'   and so should be thought of as a standard deviation parameter.
#' @inheritParams compute_taxon_specific_latent_trait_elbo
#'
#' @return A scalar. The negative ELBO associated with individual specific
#'   latent traits scaled to account for the logit transform of the heritable
#'   amplitude.
heritable_amplitude_objective <- function(
  par, i,
  heritable_amplitude, length_scale,
  taxon_specific_latent_trait_expectation,
  taxon_specific_latent_trait_outer_product_expectation,
  taxon_specific_latent_trait_covariance,
  phy,
  phylogenetic_gp,
  perform_checks = TRUE
){
  ha <- heritable_amplitude
  pgp <- phylogenetic_gp
  # tmp_par <- stats::plogis(par)
  # ha[i] <- tmp_par
  ha[i] <- par
  pgp[, , i] <- reparameterise_phylogenetic_ou(
    phy = phy,
    heritable_amplitude = par,
    length_scale = length_scale,
    environmental_amplitude = sqrt(1 - par^2),
    perform_checks = FALSE
  )
  obj <- compute_taxon_specific_latent_trait_elbo(
    taxon_specific_latent_trait_expectation = taxon_specific_latent_trait_expectation,
    taxon_specific_latent_trait_outer_product_expectation = taxon_specific_latent_trait_outer_product_expectation,
    taxon_specific_latent_trait_covariance = taxon_specific_latent_trait_covariance,
    phy = phy,
    phylogenetic_gp = pgp,
    perform_checks = perform_checks
  )
  # - (obj + log(tmp_par) + log(1 - tmp_par))
  -obj
}


#' Update ARD loading precision
#'
#' Update the ARD precision associated with each loading within the PLVM by
#' optimising the ELBO.
#'
#' @param loading_col_outer_product_expectation A D'xD'xL array of real values.
#'   The expected outer product for each column of the loading matrix under the
#'   approximate posterior distribution.
#' @param inv_loading_prior_correlation A D'xD' matrix. The inverse of the prior
#'   correlation matrix which defines the relationship between elements of the
#'   loading matrix.
#' @inheritParams initialise_loading_ard_precision
#'
#' @return An L-dimensional vector. The ARD precision associated with each
#'   column of the loading matrix.
update_loading_ard_precision <- function(
  loading_col_outer_product_expectation,
  inv_loading_prior_correlation,
  ard_shape = 1,  ard_rate = 1
){
  D <- dim(loading_col_outer_product_expectation)[1]
  L <- dim(loading_col_outer_product_expectation)[3]

  alpha <- ard_shape + (D / 2)
  beta <- sapply(
    1:L,
    function(i){
      ard_rate + (sum(diag(inv_loading_prior_correlation %*% loading_col_outer_product_expectation[, , i])) / 2)
    }
  )
  (alpha - 1) / beta
}

#' Objective function for the ordinal trait cut off point
#'
#' Update an ordinal trait cut off point by minimising the negative ELBO
#' associated with the ordinal auxiliary trait.
#'
#' @param  par A real valued scalar. The cut off point of interest
#' @param i An integer value between 3 and K. The index for the cut off point to
#'   be optimised.
#' @inheritParams compute_ordinal_auxiliary_trait_elbo
#'
#' @return  A real valued scalar. The negative ELBO associated with the ordinal
#'   auxiliary trait.
ordinal_trait_cut_off_objective <- function(
  par, i,
  y, cut_off_points,
  loading_expectation, latent_trait_expectation,
  loading_outer_expectation, latent_trait_outer_expectation,
  perform_checks = TRUE
){
  cop <- cut_off_points
  cop[i] <- par
  obj <- compute_ordinal_auxiliary_trait_elbo(
    y = y, cut_off_points = cop,
    loading_expectation = loading_expectation,
    latent_trait_expectation = latent_trait_expectation,
    loading_outer_expectation = loading_outer_expectation,
    latent_trait_outer_expectation = latent_trait_outer_expectation,
    perform_checks = perform_checks
  )
  -obj
}
