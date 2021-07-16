#' Simulate a Phylogenetic Ornstein-Uhlenbeck process
#'
#' A wrapper for \eqn{ape::rTraitCont} which simulates an Ornstein-Uhlenbeck
#' process parameterised by \eqn{amplitude > 0} and \eqn{length_scale > 0}.
#'
#' @param phy An object of class "\eqn{phylo}".
#' @inheritParams ou_kernel
#' @param internal A logical value. Specifies whether to return the values at
#'   internal nodes
#'
#' @return A numeric vector with names taken from the tip labels of \eqn{phy}. If
#'   \eqn{ancestor = TRUE}, the node labels are used if present, otherwise, “Node1”,
#'   “Node2”, etc.
#' @export
#'
#' @examples
#' tr <- ape::rtree(100)
#' z <- simulate_phylogenetic_ou(tr)
simulate_phylogenetic_ou <- function(
  phy,
  amplitude = 1, length_scale = 2,
  internal = FALSE,
  perform_checks = TRUE
  ){
  if (perform_checks) {
    checkmate::assert_number(amplitude, lower = 0, finite = TRUE)
    checkmate::assert_number(length_scale, lower = 0, finite = TRUE)
    checkmate::assert_logical(internal)
  }
  root <- stats::rnorm(1, mean = 0, sd = sqrt(amplitude))
  z <- ape::rTraitCont(
    phy = phy, model = "OU",
    sigma = sqrt(amplitude), alpha = 1 / length_scale, theta = 0,
    ancestor = internal, root.value = root
    )
  z
}

