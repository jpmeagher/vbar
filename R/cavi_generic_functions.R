#' Compute outer product expectation
#'
#' Computes the expected value of the outer product of the \eqn{D}-dimensional
#' Gaussian random variable.
#'
#' @param expected_value A real-valued vector length \eqn{L}. The expected value
#'   of the Gaussian random variable.
#' @param covariance_matrix A \eqn{L \times L} covariance matrix. The covariance
#'   matrix for the Gaussian random variable.
#' @inheritParams ou_kernel
#'
#' @return A symmetric \eqn{L \times L} matrix.
#' @export
#'
#' @examples
#' D <- 5
#' z <- rnorm(D)
#' S <- stats::rWishart(1, df = D, Sigma = diag(D)/ D^2)[, , 1]
#' gaussian_outer_product_expectation(expected_value = z, covariance_matrix = S)
gaussian_outer_product_expectation <- function(
  expected_value, covariance_matrix,
  perform_checks = TRUE
){
  if (perform_checks) {
    L <- length(expected_value)
    checkmate::assert_numeric(expected_value, any.missing = FALSE, finite = TRUE)
    checkmate::assert_matrix(covariance_matrix, mode = "numeric", any.missing = FALSE, nrows = L, ncols = L)
    checkmate::assert_true(isSymmetric(covariance_matrix))
  }
  covariance_matrix + expected_value %*% t(expected_value)
}
