#' The Ornstein-Uhlenbeck kernel
#'
#' Computes \deqn{amplitude^2 * \exp ( - d / length_scale )} where \eqn{d \ge 0}
#' is the distance between any two points in space for \eqn{amplitude > 0} and
#' \eqn{length_scale > 0}.
#'
#' @param d A set of non-negative real numbers. Pairwise distances between
#'   points.
#' @param amplitude A positive real number. The amplitude of the
#'   Ornstein-Uhlenbeck process.
#' @param length_scale A positive real number. The length-scale of the
#'   Ornstein-Uhlenbeck process.
#' @param perform_checks Logical. Check if function inputs are specified
#'   correctly.
#'
#' @return An array of positive values between 0 and \eqn{amplitude} with the
#'   same dimension as d.
#' @export
#'
#' @examples
#' x <- 1:10
#' d <- abs(outer(x, x, "-"))
#' K <- ou_kernel(d, amplitude = 1, length_scale = 2)
ou_kernel <- function(d, amplitude, length_scale, perform_checks = TRUE){
  if (perform_checks) {
    checkmate::assert_numeric(d, lower = 0, finite = TRUE, any.missing = FALSE)
    checkmate::assert_number(amplitude, lower = 0, finite = TRUE)
    checkmate::assert_number(length_scale, lower = 0, finite = TRUE)
  }
  amplitude^2 * exp(- d / length_scale)
}

#' The Exponentiated Quadratic kernel
#'
#' Computes \deqn{amplitude^2 * \exp ( - d^2 / 2 * length_scale )} where \eqn{d \ge 0}
#' is the distance between any two points in space for \eqn{amplitude > 0} and
#' \eqn{length_scale > 0}.
#'
#' @param d A set of non-negative real numbers. Pairwise distances between
#'   points.
#' @param amplitude A positive real number. The amplitude of the
#'   Gaussian process.
#' @param length_scale A positive real number. The length-scale of the
#'   Gaussian process.
#' @param perform_checks Logical. Check if function inputs are specified
#'   correctly.
#'
#' @return An array of positive values between 0 and \eqn{amplitude} with the
#'   same dimension as d.
#' @export
#'
#' @examples
#' x <- 1:10
#' d <- abs(outer(x, x, "-"))
#' K <- exp_quad_kernel(d, amplitude = 1, length_scale = 2)
exp_quad_kernel <- function(d, amplitude, length_scale, perform_checks = TRUE){
  if (perform_checks) {
    checkmate::assert_numeric(d, lower = 0, finite = TRUE, any.missing = FALSE)
    checkmate::assert_number(amplitude, lower = 0, finite = TRUE)
    checkmate::assert_number(length_scale, lower = 0, finite = TRUE)
  }
  amplitude^2 * exp(- (d^2) / ( 2 * (length_scale^2)))
}
