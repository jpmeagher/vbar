#' Initialise the Individual-specific latent traits
#'
#' Initialises individual specific latent traits given the auxiliary traits, a
#' set of loadings, and a vector of precision parameters associated with each
#' multivariate trait.
#'
#' @inheritParams compute_individual_specific_latent_trait_expectation
#' @param auxiliary_traits An NxD' matrix of real numbers. The initial matrix of
#'   auxiliary traits.
#' @param precision A positive real valued scalar. A precision value allowing the
#'   latent variables to be initialised. Larger imply that auxiliary traits are
#'   more precisely observed.
#'
#' @return An NxL dimensional matrix of real-valued latent traits.
initialise_individual_specific_latent_traits <- function(
  auxiliary_traits,
  loading,
  precision = 1,
  perform_checks = TRUE
){
  D <- ncol(auxiliary_traits)
  N <- nrow(auxiliary_traits)
  L <- ncol(loading)
  if (perform_checks) {
    checkmate::assert_matrix(
      auxiliary_traits, mode = "numeric",
      nrows = N, ncols = D, any.missing = FALSE
    )
    checkmate::assert_matrix(
      loading, mode = "numeric",
      nrows = D, ncols = L, any.missing = FALSE
    )
  }
  M_inv <- chol2inv(chol(t(loading) %*% loading + (1 / precision) * diag(L)))
  Z <- auxiliary_traits %*% (loading %*% M_inv)
  Z
}



#' Compute Individual-specific Latent Trait Precision
#'
#' Computes the precision matrix for the approximate posterior distribution of
#' individual-specific latent traits within the PLVM.
#'
#' @param precision_vector A D-dimensional vector of positive
#'   real values The precision parameters associated with each of the D
#'   auxiliary traits.
#' @param loading_outer_product_expectation An \eqn{LxLxD} dimensional array of
#'   real values. The expected outer product of the loading under the
#'   approximate posterior distribution.
#' @param within_taxon_amplitude An L-dimensional vector of real values
#'   on the interval \eqn{[0, 1]}. The within-taxon variation parameters for
#'   each of the L latent traits.
#' @inheritParams ou_kernel
#'
#' @return A symmetric LxL matrix of real values. The precision matrix of the
#'   individual-specific latent traits under the approximate posterior
#'   distribution.
compute_individual_specific_latent_trait_precision <- function(
  precision_vector,
  loading_outer_product_expectation,
  within_taxon_amplitude,
  perform_checks = TRUE
){
  D <- length(precision_vector)
  L <- length(within_taxon_amplitude)
  if (perform_checks) {
    checkmate::assert_numeric(
      precision_vector, any.missing = FALSE
    )
    checkmate::assert_array(
      loading_outer_product_expectation, mode = "numeric", any.missing = FALSE,
      d = 3
    )
    checkmate::assert_true(
      all(dim(loading_outer_product_expectation) == c(L, L, D))
    )
    checkmate::assert_numeric(
      within_taxon_amplitude, lower = 0, upper = 1, any.missing = FALSE
    )
  }
  apply(
    sweep(loading_outer_product_expectation, 3, precision_vector, "*"),
    c(1, 2), sum
    ) + (diag(L) *(1 / within_taxon_amplitude^2))
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
compute_individual_specific_latent_trait_expectation <- function(
  auxiliary_trait,
  loading,
  taxon_specific_latent_trait,
  precision_vector,
  individual_specific_latent_trait_precision,
  within_taxon_amplitude,
  perform_checks = TRUE
){
  D <- length(precision_vector)
  L <- length(within_taxon_amplitude)
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
      precision_vector, any.missing = FALSE
    )
    checkmate::assert_matrix(
      individual_specific_latent_trait_precision, mode = "numeric",
      nrows = L, ncols = L, any.missing = FALSE
    )
    checkmate::assert_numeric(
      within_taxon_amplitude, lower = 0, upper = 1, any.missing = FALSE
    )
  }
  tmp <- t(loading) %*% diag(precision_vector) %*% auxiliary_trait +
    (diag(L)*(1 / within_taxon_amplitude^2)) %*% taxon_specific_latent_trait
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
#' @param N_s A natural number. The number of individuals within the taxon.
#' @inheritParams  compute_individual_specific_latent_trait_precision
#' @param conditional_standard_deviation An L-dimensional vector of real values
#'   on the interval \eqn{[0, 1]}. The standard deviation of the taxon-specific
#'   latent trait conditional on its parent taxon.
#'
#' @return An L-dimensional vector of real values. Precision parameters for the
#'   taxo-specific latent traits.
compute_terminal_taxon_specific_latent_trait_precision <- function(
  N_s,
  within_taxon_amplitude,
  conditional_standard_deviation,
  perform_checks = TRUE
){
  L <- length(within_taxon_amplitude)
  if (perform_checks) {
    checkmate::assert_numeric(
      within_taxon_amplitude, lower = 0, upper = 1, any.missing = FALSE
    )
    checkmate::assert_numeric(
      conditional_standard_deviation, lower = 0, upper = 1, any.missing = FALSE
    )
  }
  N_s / within_taxon_amplitude^2 +
    1 / conditional_standard_deviation^2
}

#' Compute Terminal Taxon-specific Latent Trait Expectation
#'
#' Computes the expected value of taxon-specific latent traits at terminal nodes
#' of the phylogeny under the approximate posterior distribution for the PLVM.
#'
#' @param individual_specific_latent_traits An L-dimensional vector of real
#'   values. The latent trait values associated with individuals within the
#'   taxon.
#' @inheritParams compute_terminal_taxon_specific_latent_trait_precision
#' @param parent_taxon_latent_trait An L-dimensional vector of real values. The
#'   latent trait values associated with the parent taxon.
#' @param conditional_expectation_weight An L-dimensional vector of real values
#'   on the interval \eqn{[0, 1]}. The weight of the latent trait values for the
#'   parent taxon in the conditional distribution for the taxon-specific latent
#'   trait.
#' @param latent_trait_precision An L-dimensional vector of positive real
#'   values. The precision of the taxon-specific latent trait.
#'
#' @return An L-dimensional vector of real values.
compute_terminal_taxon_specific_latent_trait_expectation <- function(
  individual_specific_latent_traits,
  within_taxon_amplitude,
  parent_taxon_latent_trait,
  conditional_expectation_weight,
  conditional_standard_deviation,
  latent_trait_precision,
  perform_checks = TRUE
){
  L <- length(within_taxon_amplitude)
  if (perform_checks) {
    checkmate::assert_matrix(
      individual_specific_latent_traits, mode = "numeric",
      any.missing = FALSE, ncols = L
    )
    checkmate::assert_numeric(
      within_taxon_amplitude, lower = 0, upper = 1, any.missing = FALSE
    )
    checkmate::assert_numeric(
      parent_taxon_latent_trait, any.missing = FALSE, len = L
    )
    checkmate::assert_numeric(
      conditional_expectation_weight, lower = 0, upper = 1, any.missing = FALSE, len = L
    )
    checkmate::assert_numeric(
      conditional_standard_deviation, lower = 0, upper = 1, any.missing = FALSE, len = L
    )
    checkmate::assert_numeric(
      latent_trait_precision, any.missing = FALSE, len = L, lower = 0
    )
  }
  tmp <- colSums(sweep(individual_specific_latent_traits, 2, within_taxon_amplitude^2, "/")) +
    conditional_expectation_weight * parent_taxon_latent_trait / conditional_standard_deviation^2
  tmp / latent_trait_precision
}

#' Compute Internal Taxon-specific Latent Trait Precision
#'
#' Computes the precision of taxon-specific latent traits at internal nodes of
#' the phylogeny under the approximate posterior distribution for the PLVM.
#'
#' @param child_taxa_conditional_expectation_weights An L-dimensional vector of
#'   real values on the interval \eqn{[0, 1]}. The weight of the taxon-specific latent trait
#'   values in the conditional distribution of the latent
#'   traits for its child taxa.
#' @param child_taxa_conditional_standard_deviations An L-dimensional vector of
#'   real values on the interval \eqn{[0, 1]}. The standard deviation of the taxon-specific
#'   latent traits of child taxa conditional on the parent taxon.
#' @inheritParams compute_terminal_taxon_specific_latent_trait_precision
#'
#' @return An L-dimensional vector of positive real values.
compute_internal_taxon_specific_latent_trait_precision <- function(
  child_taxa_conditional_expectation_weights,
  child_taxa_conditional_standard_deviations,
  conditional_standard_deviation,
  perform_checks = TRUE
){
  L <- length(conditional_standard_deviation)
  if (perform_checks) {
    checkmate::assert_matrix(
      child_taxa_conditional_expectation_weights, mode = "numeric",
      any.missing = FALSE, ncols = L
    )
    checkmate::assert_numeric(
      child_taxa_conditional_expectation_weights,
      lower = 0, upper = 1
    )
    checkmate::assert_matrix(
      child_taxa_conditional_standard_deviations, mode = "numeric",
      any.missing = FALSE, ncols = L
    )
    checkmate::assert_numeric(
      child_taxa_conditional_standard_deviations,
      lower = 0, upper = 1
    )
    checkmate::assert_numeric(
      conditional_standard_deviation, lower = 0, upper = 1, any.missing = FALSE, len = L
    )
  }
  colSums(child_taxa_conditional_expectation_weights^2 / child_taxa_conditional_standard_deviations^2) +
    1 / conditional_standard_deviation^2
}

#' Compute Internal Taxon-specific Latent Trait Expectation
#'
#' Computes the expected value of taxon-specific latent traits at internal nodes
#' of the phylogeny under the approximate posterior distribution for the PLVM.
#'
#' @param child_taxa_latent_traits An L-dimensional vector of real values. The latent trait values associated with child taxa.
#' @inheritParams compute_internal_taxon_specific_latent_trait_precision
#' @inheritParams compute_terminal_taxon_specific_latent_trait_expectation
#'
#' @return An L-dimensional vector of real values.
compute_internal_taxon_specific_latent_trait_expectation <- function(
  child_taxa_latent_traits,
  child_taxa_conditional_expectation_weights,
  child_taxa_conditional_standard_deviations,
  parent_taxon_latent_trait,
  conditional_expectation_weight,
  conditional_standard_deviation,
  latent_trait_precision,
  perform_checks = TRUE
){
  L <- length(conditional_standard_deviation)
  if (perform_checks) {
    checkmate::assert_matrix(
      child_taxa_latent_traits, mode = "numeric",
      any.missing = FALSE, ncols = L
    )
    checkmate::assert_matrix(
      child_taxa_conditional_expectation_weights, mode = "numeric",
      any.missing = FALSE, ncols = L
    )
    checkmate::assert_numeric(
      child_taxa_conditional_expectation_weights,
      lower = 0, upper = 1
    )
    checkmate::assert_matrix(
      child_taxa_conditional_standard_deviations, mode = "numeric",
      any.missing = FALSE, ncols = L
    )
    checkmate::assert_numeric(
      child_taxa_conditional_standard_deviations,
      lower = 0, upper = 1
    )
    checkmate::assert_numeric(
      parent_taxon_latent_trait, any.missing = FALSE, len = L
    )
    checkmate::assert_numeric(
      conditional_expectation_weight, lower = 0, upper = 1, any.missing = FALSE, len = L
    )
    checkmate::assert_numeric(
      conditional_standard_deviation, lower = 0, upper = 1, any.missing = FALSE, len = L
    )
    checkmate::assert_numeric(
      latent_trait_precision, any.missing = FALSE, len = L, lower = 0
    )
  }
  tmp <- colSums(child_taxa_conditional_expectation_weights * child_taxa_latent_traits / child_taxa_conditional_standard_deviations^2) +
    conditional_expectation_weight * parent_taxon_latent_trait / conditional_standard_deviation^2
  tmp / latent_trait_precision
}

#' Individual Specific Latent Trait Trait ELBO
#'
#' Compute the contribution of individual specific latent traits to the Evidence
#' Lower Bound (ELBO) of a PLVM given the approximate posterior distribution for
#' loadings and latent traits
#'
#' @param individual_specific_latent_trait_expectation A NxL matrix of real
#'   values. The expected value of the individual specific latent traits.
#' @param taxon_id A N-dimensional vector of labels. Matches each individual
#'   specific latent trait to a terminal node in \eqn{phy}. Each element must
#'   correspond to one of \eqn{phy$tip.labels}.
#' @inheritParams simulate_phylogenetic_ou
#' @param terminal_taxon_specific_latent_trait_expectation A SxL matrix of real
#'   values. The expected taxon-specific latent trait for terminal nodes in
#'   \eqn{phy}. Rows must be ordered by \eqn{phy$tip.labels}.
#' @param individual_specific_latent_trait_covariance A LxL covariance matrix.
#'   The covariance of individual specific latent traits under the approximate
#'   posterior.
#' @param individual_specific_latent_trait_outer_product_expectation A LxLxN
#'   array of real values. The expected outer product of individual specific
#'   latent traits under the approximate posterior.
#' @param terminal_taxon_latent_trait_outer_product_expectation A LxLxS array of
#'   real values. The expected outer product of taxon specific latent traits at
#'   terminal nodes under the approximate posterior.
#' @inheritParams compute_individual_specific_latent_trait_precision
#'
#' @return A real valued scalar. The contribution of individual specific latent
#'   traits to the ELBO.
compute_individual_specific_latent_trait_elbo <- function(
  individual_specific_latent_trait_expectation,
  taxon_id, phy,
  terminal_taxon_specific_latent_trait_expectation,
  individual_specific_latent_trait_covariance,
  individual_specific_latent_trait_outer_product_expectation,
  terminal_taxon_latent_trait_outer_product_expectation,
  within_taxon_amplitude,
  perform_checks = TRUE
){
  N <- nrow(individual_specific_latent_trait_expectation)
  L <- ncol(individual_specific_latent_trait_expectation)
  S <- length(phy$tip.label)
  if (perform_checks) {
    checkmate::assert_matrix(individual_specific_latent_trait_expectation, mode = "numeric", any.missing = FALSE)
    checkmate::assert_set_equal(
      sort(phy$tip.label),
      sort(taxon_id)
    )
    checkmate::assert_matrix(
      terminal_taxon_specific_latent_trait_expectation, mode = "numeric", any.missing = FALSE,
      nrows = S, ncol = L
    )
    checkmate::assert_matrix(
      individual_specific_latent_trait_covariance, mode = "numeric", any.missing = FALSE,
      nrows = L, ncol = L
    )
    checkmate::assert_array(
      terminal_taxon_latent_trait_outer_product_expectation, mode = "numeric",
      any.missing = FALSE, d = 3
    )
    checkmate::assert_set_equal(
      dim(terminal_taxon_latent_trait_outer_product_expectation),
      c(L, L, S)
    )
    checkmate::assert_array(
      individual_specific_latent_trait_outer_product_expectation, mode = "numeric",
      any.missing = FALSE, d = 3
    )
    checkmate::assert_set_equal(
      dim(individual_specific_latent_trait_outer_product_expectation),
      c(L, L, N)
    )
    checkmate::assert_numeric(
      within_taxon_amplitude, len = L, any.missing = FALSE
    )
  }
  individuals_per_taxon <- as.numeric(
    table(taxon_id)[phy$tip.label]
    )
  A1 <- - N * sum(log(within_taxon_amplitude))
  A2.1 <- apply(
    individual_specific_latent_trait_outer_product_expectation,
    c(1, 2), sum
  )
  A2.2 <- apply(sweep(
    terminal_taxon_latent_trait_outer_product_expectation,
    3, individuals_per_taxon, "*"
    ), c(1, 2), sum
  )
  A2 <- -0.5 * sum(diag(
    (A2.1 + A2.2) %*% (diag(L) * (within_taxon_amplitude^-2))
  ))
  A3 <- sum(sapply(
    1:S,
    function(i){
      t(terminal_taxon_specific_latent_trait_expectation[i, ]) %*%
        (diag(L) * (within_taxon_amplitude^-2)) %*%
        colSums(
          individual_specific_latent_trait_expectation[taxon_id == phy$tip.label[i], ,drop=F]
        )
    }
  ))
  A4 <- N * sum(log(diag(chol(individual_specific_latent_trait_covariance))))
  A1 + A2 + A3 + A4
}


#' Taxon Specific Latent Trait Trait ELBO
#'
#' Compute the contribution of taxon specific latent traits to the Evidence
#' Lower Bound (ELBO) of a PLVM given the approximate posterior distribution for
#' loadings and latent traits
#'
#' @inheritParams simulate_phylogenetic_ou
#' @param taxon_specific_latent_trait_expectation A (2S - 1)xL matrix of real
#'   values. The expected value of the taxon specific latent traits.
#' @param taxon_specific_latent_trait_outer_product_expectation A LxLx(2S - 1)
#'   array of real values. The expected outer product of taxon specific latent
#'   traits under the approximate posterior.
#' @param taxon_specific_latent_trait_covariance A LxLx(2S - 1) array of real
#'   values. The covariance of taxon specific latent traits under the
#'   approximate posterior.
#' @param phylogenetic_gp A (2S - 1)x2xL array of real values. Hyperparameters
#'   for the Gaussian process over \eqn{phy}.
#'
#' @return A real valued scalar. The contribution of taxon specific latent
#'   traits to the ELBO.
compute_taxon_specific_latent_trait_elbo <- function(
  taxon_specific_latent_trait_expectation,
  taxon_specific_latent_trait_outer_product_expectation,
  taxon_specific_latent_trait_covariance,
  phy,
  phylogenetic_gp,
  perform_checks = TRUE
){
  S <- length(phy$tip.label)
  L <- ncol(taxon_specific_latent_trait_expectation)
  if (perform_checks) {
    checkmate::assert_matrix(
      taxon_specific_latent_trait_expectation, mode = "numeric", any.missing = FALSE,
      nrow = 2*S - 1
    )
    checkmate::assert_array(
      taxon_specific_latent_trait_outer_product_expectation,
      mode = "numeric", any.missing = FALSE, d = 3
    )
    checkmate::assert_set_equal(
      dim(taxon_specific_latent_trait_outer_product_expectation),
      c(L, L, 2*S - 1)
    )
    checkmate::assert_matrix(taxon_specific_latent_trait_covariance)
    checkmate::assert_array(
      phylogenetic_gp,
      mode = "numeric", any.missing = FALSE, d = 3
    )
    checkmate::assert_set_equal(
      dim(phylogenetic_gp),
      c(2*S - 1, 2, L)
    )
  }
  traversal_order <- phy$edge[ape::postorder(phy), ]

  A1 <- - sum(log(phylogenetic_gp[, "sd", ]))
  A2.1 <- sum(sapply(
    1:(2 *S - 2),
    function(i) {
      sum(diag(
        (taxon_specific_latent_trait_outer_product_expectation[, , traversal_order[i, 2]] +
          ((diag(L) * (phylogenetic_gp[traversal_order[i, 2], "weight", ]^2)) %*%
             taxon_specific_latent_trait_outer_product_expectation[, , traversal_order[i, 1]])
        ) %*%
          (diag(L) * (phylogenetic_gp[traversal_order[i, 2], "sd", ]^-2))
      ))
    }
  ))
  A2.2 <- sum(diag(
    taxon_specific_latent_trait_outer_product_expectation[, , S + 1] %*%
      (diag(L) * (phylogenetic_gp[S + 1, "sd", ]^-2))
  ))
  A2 <- - 0.5 * (A2.1 + A2.2)
  A3 <- sum(sapply(
    1:(2 * S - 2),
    function(i){
      t(taxon_specific_latent_trait_expectation[traversal_order[i, 1], ] *
          phylogenetic_gp[traversal_order[i, 2], "weight", ] *
          (phylogenetic_gp[traversal_order[i, 2], "sd", ]^-2)) %*%
        taxon_specific_latent_trait_expectation[traversal_order[i, 2], ]
    }
    ))
  A4 <- 0.5 * sum(log(exp(1) * taxon_specific_latent_trait_covariance))
  A1 + A2 + A3 + A4
}


compute_latent_trait_elbo <- function(
  individual_specific_latent_trait_expectation,
  taxon_id, phy,
  terminal_taxon_specific_latent_trait_expectation,
  individual_specific_latent_trait_covariance,
  individual_specific_latent_trait_outer_product_expectation,
  terminal_taxon_latent_trait_outer_product_expectation,
  taxon_specific_latent_trait_expectation,
  taxon_specific_latent_trait_outer_product_expectation,
  taxon_specific_latent_trait_covariance,
  phylogenetic_gp,
  within_taxon_amplitude,
  perform_checks = TRUE
){
  elbo <- compute_individual_specific_latent_trait_elbo(
    individual_specific_latent_trait_expectation = individual_specific_latent_trait_expectation,
    taxon_id = taxon_id, phy = phy,
    terminal_taxon_specific_latent_trait_expectation = terminal_taxon_specific_latent_trait_expectation,
    individual_specific_latent_trait_covariance = individual_specific_latent_trait_covariance,
    individual_specific_latent_trait_outer_product_expectation = individual_specific_latent_trait_outer_product_expectation,
    terminal_taxon_latent_trait_outer_product_expectation = terminal_taxon_latent_trait_outer_product_expectation,
    within_taxon_amplitude = within_taxon_amplitude,
    perform_checks = perform_checks
  ) +
    compute_taxon_specific_latent_trait_elbo(
      taxon_specific_latent_trait_expectation = taxon_specific_latent_trait_expectation,
      taxon_specific_latent_trait_outer_product_expectation = taxon_specific_latent_trait_outer_product_expectation,
      taxon_specific_latent_trait_covariance = taxon_specific_latent_trait_covariance,
      phy = phy,
      phylogenetic_gp = phylogenetic_gp,
      perform_checks = perform_checks
    )
  elbo
}
