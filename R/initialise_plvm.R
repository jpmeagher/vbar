#' Initialise CAVI for a PLVM
#'
#' Specify a set of initial values for Co-ordinate Ascent Variational Inference
#' (CAVI) on the Phylogenetic Latent Variable Model (PLVM).
#'
#' @inheritParams initialise_auxiliary_traits
#' @inheritParams initialise_loading
#' @inheritParams initialise_loading_ard_precision
#' @inheritParams initialise_precision
#' @inheritParams reparameterise_phylogenetic_ou
#' @inheritParams compute_individual_specific_latent_trait_precision
#' @param metadata A data frame. Contains all the metadata required to map a set
#'   of manifest traits to the PLVM auxiliary traits..
#' @param id_label Either a string or numeric indexing the variable in
#'   \eqn{manifest_trait_df} identifying the taxon to which each observation
#'   belongs.
#' @param random_seed A single value, interpreted as an integer, or
#'   NULL.
#'
#' @seealso specify_manifest_trait_metadata
#'
#' @return A list of intital values for a PLVM..
#' @export
initialise_plvm <- function(
  manifest_trait_df, metadata, L, phy, id_label = "taxon_id",
  loading_prior_correlation,
  auxiliary_traits = NULL,
  precision = NULL,
  ard_precision = rep(1, L),
  ard_shape = 1, ard_rate = 1,
  loading = NULL, method = "varimax",
  within_taxon_amplitude = NULL,
  heritable_amplitude = NULL,
  length_scale = 2,
  perform_checks = TRUE,
  random_seed = NULL
){
  # Check Metadata and manifest traits
  metadata <- specify_manifest_trait_metadata(
    n_traits = nrow(metadata),
    trait_names = metadata$trait_names,
    trait_type = metadata$trait_type,
    trait_levels = metadata$trait_levels,
    manifest_trait_index = metadata$manifest_trait_index,
    auxiliary_trait_index = metadata$auxiliary_trait_index,
    link_functions = metadata$link_functions,
    inverse_link_functions = metadata$inverse_link_functions,
    cut_off_points = metadata$cut_off_points,
    categories = metadata$categories,
    manifest_trait_df = manifest_trait_df,
    perform_checks = perform_checks
  )
  if (perform_checks) {
    checkmate::assert_set_equal(
      sort(phy$tip.label),
      sort(unique(manifest_trait_df[, id_label]))
    )
  }
  set.seed(random_seed)
  # Set Dimensions
  P <- nrow(metadata)
  N <- nrow(manifest_trait_df)
  S <- length(phy$tip.label)
  D <- sum(sapply(metadata$manifest_trait_index, length))
  D_prime <- sum(sapply(metadata$auxiliary_trait_index, length))
  # Hyper-parameters
  alpha <- initialise_loading_ard_precision(
    L = L,
    ard_shape = ard_shape, ard_rate = ard_rate,
    ard_precision = ard_precision,
    perform_checks = perform_checks
  )
  if (is.null(within_taxon_amplitude)) within_taxon_amplitude <- 0.05 + 0.1 * stats::runif(L)
  if (is.null(heritable_amplitude)) heritable_amplitude <- 0.5 + 0.25 * stats::runif(L)
  if (perform_checks) {
    checkmate::assert_numeric(
      heritable_amplitude, lower = 0, upper = 1, any.missing = FALSE, len = L
    )
  }
  phylogenetic_GP <- lapply(
    1:L, function(i){
      reparameterise_phylogenetic_ou(
        phy = phy,
        heritable_amplitude = heritable_amplitude[i],
        length_scale = length_scale,
        environmental_amplitude = sqrt(1 - heritable_amplitude[i]^2),
        perform_checks = perform_checks
      )
    }
  )
  phylogenetic_GP <- simplify2array(phylogenetic_GP)
  # Auxiliary Traits
  X <- initialise_auxiliary_traits(
    n_traits = nrow(metadata),
    manifest_trait_df = manifest_trait_df,
    trait_names = metadata$trait_names,
    trait_type = metadata$trait_type,
    trait_levels = metadata$trait_levels,
    manifest_trait_index = metadata$manifest_trait_index,
    auxiliary_trait_index = metadata$auxiliary_trait_index,
    inverse_link_functions = metadata$inverse_link_functions,
    auxiliary_traits = auxiliary_traits,
    perform_checks = perform_checks
  )
  # Latent Variables
  W <- initialise_loading(
    D_prime = D_prime, L = L,
    ard_precision = alpha,
    loading_prior_correlation = loading_prior_correlation,
    loading = loading, method = method, random_seed = random_seed,
    auxiliary_traits = X,
    perform_checks = perform_checks
  )
  Z <- initialise_individual_specific_latent_traits(
    auxiliary_traits = X,
    loading = W, precision = 1,
    perform_checks = perform_checks
  )
  # Precision Parameters
  lambda <- initialise_precision(
    n_traits = P,
    trait_names = metadata$trait_names,
    trait_type = metadata$trait_type,
    precision_prior_shape = 100,
    precision_prior_rate = 100,
    precision = precision,
    perform_checks = perform_checks
  )
  lambda_vector <- map_precision_to_auxiliary_traits(
    precision = lambda,
    auxiliary_trait_index = metadata$auxiliary_trait_index,
    perform_checks = perform_checks
  )
  # Loading Parameters
  c_star <- compute_scaled_conditional_row_variance_vector(
    loading_prior_correlation
    )
  loading_row_conditional_mean_weight <- compute_loading_row_conditional_mean_weight_matrix(
    loading_prior_correlation
    )
  U_C_w <- chol(loading_prior_correlation)
  inv_C_w <- chol2inv(U_C_w)
  log_det_C_w <- 2 * sum(log(diag(U_C_w)))
  lambda_W_list <- compute_loading_row_precision_list(
    total_individual_specific_latent_trait_outer_product_expectation = t(Z) %*% Z,
    precision_vector = lambda_vector,
    ard_precision = alpha,
    scaled_conditional_row_variance_vector = c_star,
    perform_checks = perform_checks
  )
  lambda_W <- simplify2array(lambda_W_list)
  lambda_W <- array(lambda_W, dim = c(L, L, D_prime))
  inv_lambda_W <- simplify2array(lapply(
    1:D_prime,
    function(i){
      chol2inv(chol(lambda_W[, , i]))
    }
  ))
  inv_lambda_W <- array(inv_lambda_W, dim = c(L, L, D_prime))
  outer_W_list <- lapply(
    1:D_prime,
    function(i){
      gaussian_outer_product_expectation(
        expected_value = W[i, ], covariance_matrix = matrix(inv_lambda_W[, , i], nrow = L, ncol = L),
        perform_checks = perform_checks
      )
    }
  )
  outer_W <- simplify2array(outer_W_list)
  outer_W <- array(outer_W, dim = c(L, L, D_prime))
  outer_W_col_list <- lapply(
    1:L,
    function(i){
      gaussian_outer_product_expectation(
        expected_value = W[, i], covariance_matrix = diag(inv_lambda_W[i, i, ]),
        perform_checks = perform_checks
      )
    }
  )
  outer_W_col <- simplify2array(outer_W_col_list)
  # Individual Specific Latent Trait Parameters
  tmp_W_outer <- simplify2array(outer_W_list)
  tmp_W_outer <- array(tmp_W_outer, dim = c(L, L, D_prime))
  lambda_Z <- compute_individual_specific_latent_trait_precision(
    precision_vector = lambda_vector,
    loading_outer_product_expectation = tmp_W_outer,
    within_taxon_amplitude = within_taxon_amplitude,
    perform_checks = perform_checks
  )
  inv_lambda_Z <- chol2inv(chol(lambda_Z))
  outer_Z_list <- lapply(
    1:N,
    function(i){
      gaussian_outer_product_expectation(
        expected_value = Z[i, ], covariance_matrix = inv_lambda_Z,
        perform_checks = perform_checks
      )
    }
  )
  outer_Z <- simplify2array(outer_Z_list)
  outer_Z <- array(outer_Z, dim = c(L, L, N))
  # Taxon Specific Latent Traits
  f <- lambda_F <- matrix(NA, S + phy$Nnode, L)
  outer_F <- array(NA, dim = c(L, L, S + phy$Nnode))
  # Terminal nodes
  for (i in 1:S) {
    lambda_F[i, ] <- compute_terminal_taxon_specific_latent_trait_precision(
      N_s = sum(manifest_trait_df[, id_label] == phy$tip.label[i]),
      within_taxon_amplitude = within_taxon_amplitude,
      conditional_standard_deviation = phylogenetic_GP[i, "sd", ],
      perform_checks = perform_checks
    )
    f[i, ] <- colMeans(as.matrix(Z[manifest_trait_df[, id_label] == phy$tip.label[i], ]))
    outer_F[, , i] <- gaussian_outer_product_expectation(
      expected_value = f[i, ], precision_matrix = diag(L)  * lambda_F[i, ],
      perform_checks = perform_checks
    )
  }
  for (i in unique(phy$edge[ape::postorder(phy), 1])) {
    ch <- phy$edge[phy$edge[, 1] == i, 2]
    lambda_F[i, ] <- compute_internal_taxon_specific_latent_trait_precision(
      child_taxa_conditional_expectation_weights = as.matrix(phylogenetic_GP[ch, "weight", ]),
      child_taxa_conditional_standard_deviations = as.matrix(phylogenetic_GP[ch, "sd", ]),
      conditional_standard_deviation = phylogenetic_GP[i, "sd", ],
      perform_checks = perform_checks
    )
    f[i, ] <- colMeans(as.matrix(f[ch, ]))
    outer_F[, , i] <- gaussian_outer_product_expectation(
      expected_value = f[i, ], covariance_matrix = diag(L) * (1 / lambda_F[i, ]),
      perform_checks = perform_checks
    )
  }

  out <- list(
    manifest_trait_df = manifest_trait_df,
    metadata = metadata,
    phy = phy, id_label = id_label,
    auxiliary_traits = X,
    precision = lambda,
    precision_vector = lambda_vector,
    ard_precision = alpha,
    scaled_conditional_loading_row_variance_vector = c_star,
    loading_row_conditional_mean_weight = loading_row_conditional_mean_weight,
    loading_prior_correlation = loading_prior_correlation,
    inv_loading_prior_correlation = inv_C_w,
    loading_prior_correlation_log_det = log_det_C_w,
    loading_expectation = W,
    loading_row_precision = lambda_W,
    loading_row_covariance = inv_lambda_W,
    loading_row_outer_product_expectation = outer_W,
    loading_col_outer_product_expectation = outer_W_col,
    within_taxon_amplitude = within_taxon_amplitude,
    individual_specific_latent_trait_precision = lambda_Z,
    individual_specific_latent_trait_covariance = inv_lambda_Z,
    individual_specific_latent_trait_expectation = Z,
    individual_specific_latent_trait_outer_product_expectation = outer_Z,
    heritable_amplitude = heritable_amplitude,
    length_scale = length_scale,
    phylogenetic_GP = phylogenetic_GP,
    taxon_specific_latent_trait_precision = lambda_F,
    taxon_specific_latent_trait_covariance = 1 / lambda_F,
    taxon_specific_latent_trait_expectation = f,
    taxon_specific_latent_trait_outer_product_expectation = outer_F
  )

  out
}
