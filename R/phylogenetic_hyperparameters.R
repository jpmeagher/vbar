#' Phylogenetic Ornstein-Uhlenbeck hyperparameters
#'
#' Compute the hyper-parameters at each node on a phylogeny given the amplitude
#' and length-scale of the Ornstein-Uhlenbeck process.
#'
#' @inheritParams simulate_phylogenetic_ou
#'
#' @return A \eqn{2n - 1 \times 2} matrix of positive real values, where \eqn{n}
#'   the the number of terminal nodes on a bifurcating tree. The first column is
#'   the contribution of the parent node to the expected value, while the second
#'   is the conditional variance given the parent.
#' @export
#'
#' @examples
#' tr <- ape::rtree(100)
#' hyper <- phylogenetic_ou_hyperparameters(tr)
phylogenetic_ou_hyperparameters <- function(
  phy,
  amplitude = 1, length_scale = 2,
  perform_checks = TRUE
){
  if (perform_checks) {
    checkmate::assert_number(amplitude, lower = 0, finite = TRUE)
    checkmate::assert_number(length_scale, lower = 0, finite = TRUE)
    if (is.null(phy$edge.length))
      stop("tree has no branch length")
    if (any(phy$edge.length < 0))
      stop("at least one branch length negative")
  }
  n <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- n + 1L
  hyper <- matrix(nrow = n + phy$Nnode, ncol = 2)
  hyper[ROOT, ] <- c(1, amplitude^2)
  des <- phy$edge[, 2]
  d <- phy$edge.length
  for (i in 1:N) {
    hyper[des[i], ] <- c(exp(- d[i] / length_scale), amplitude^2 * (1 - exp(- 2 * d[i] / length_scale)) )
  }
  colnames(hyper) <- c("parent_weight", "variance")
  hyper
}
