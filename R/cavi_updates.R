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
          sum(diag(2 * t(X[, ind]) %*% M[, ind])) +
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

update_phylogenetic_gp <- function(){}
