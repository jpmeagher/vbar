#' Compute latent trait outer product expectations
#'
#' Computes the expected value of the outer product of the \eqn{L} latent traits for the \eqn{n^{th}} individual under the approximate posterior distribution.
#'
#' @param z_n_tilde A real-valued vector length \eqn{L}. The expected latent trait under the approximate posterior distribution.
#' @param S_n_z A \eqn{L \times L} covariance matrix. The variance of the latent trait under the approximate posterior distribution.
#' @inheritParams ou_kernel
#'
#' @return A symmetric \eqn{L \times L} matrix.
#' @export
#'
#' @examples
#' L <- 5
#' z <- rnorm(L)
#' S <- stats::rWishart(1, df = L, Sigma = diag(L)/ L^2)[, , 1]
#' latent_trait_outer_product_expectation(z_n_tilde = z, S_n_z = S)
latent_trait_outer_product_expectation <- function(
  z_n_tilde, S_n_z,
  perform_checks = TRUE
){
  if (perform_checks) {
    L <- length(z_n_tilde)
    checkmate::assert_numeric(z_n_tilde, any.missing = FALSE, finite = TRUE)
    checkmate::assert_matrix(S_n_z, mode = "numeric", any.missing = FALSE, nrows = L, ncols = L)
    checkmate::assert_true(isSymmetric(S_n_z))
  }
  S_n_z + z_n_tilde %*% t(z_n_tilde)
}

