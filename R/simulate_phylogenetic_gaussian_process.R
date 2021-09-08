#' Simulate a Phylogenetic Ornstein-Uhlenbeck process
#'
#' A wrapper for \eqn{ape::rTraitCont} which simulates an Ornstein-Uhlenbeck
#' process over the phylogeny \eqn{phy} such that the phylogenetic Gaussian
#' process \deqn{ f ( \boldsymbol t )) \sim \mathcal{GP} ( 0, k (t, t' ))))} for
#' \deqn{ k (t, t' )) = heritable_amplitude^2 \exp ( - \frac{d (t, t'
#' ))}{length_scale})) + environmental_amplitude^2 \delta ( t_1, \dots, t_S )) }
#' parameterised by \eqn{heritable_amplitude > 0}, \eqn{length_scale > 0}, and
#' \eqn{environmental_amplitude > 0}, where \eqn{t} is a position on the
#' phylogeny, \eqn{d (t, t')} is the length of the shortest path from \eqn{t} to
#' \eqn{t'}, \eqn{t_1, \dots, t_S} are positions corresponding to the \eqn{S}
#' extant taxa, and \eqn{\delta ( \cdot ))} is the indicator function.
#'
#' @param phy An object of class "\eqn{phylo}".
#' @param heritable_amplitude A positive real-valued scalar. The amplitude of
#'   the phylogenetic OU process. Note that this is the square root of the
#'   variance and so should be thought of as a standard deviation parameter.
#' @param environmental_amplitude A positive real-valued scalar. The amplitude
#'   of non-phylogenetic variation in the traits of extant taxa. N Note that
#'   this is the square root of the variance and so should be thought of as a
#'   standard deviation parameter.
#' @inheritParams ou_kernel
#' @param internal A logical value. Specifies whether to return the values at
#'   internal nodes
#'
#' @return A numeric vector with names taken from the tip labels of \eqn{phy}.
#'   If \eqn{ancestor = TRUE}, the node labels are used if present, otherwise,
#'   “Node1”, “Node2”, etc.
#' @export
#'
#' @examples
#' S <- 100
#' tr <- ape::rtree(S)
#' z <- simulate_phylogenetic_ou(tr)
simulate_phylogenetic_ou <- function(
  phy,
  heritable_amplitude = 1, length_scale = 2,
  environmental_amplitude = 0,
  internal = FALSE,
  perform_checks = TRUE
  ){
  if (perform_checks) {
    checkmate::assert_number(heritable_amplitude, lower = 0, finite = TRUE)
    checkmate::assert_number(length_scale, lower = 0, finite = TRUE)
    checkmate::assert_number(environmental_amplitude, lower = 0, finite = TRUE)
    checkmate::assert_logical(internal)
  }
  S <- length(phy$tip.label)
  root <- stats::rnorm(1, mean = 0, sd = heritable_amplitude)
  f <- ape::rTraitCont(
    phy = phy, model = "OU",
    sigma = heritable_amplitude, alpha = 1 / length_scale, theta = 0,
    ancestor = TRUE, root.value = root
    )
  f[1:S] <- f[1:S] + stats::rnorm(S, sd = environmental_amplitude)
  if (internal) {
    return(f)
  } else {
    return(f[1:S])
  }
}

