#' Reparameterise phylogenetic Ornstein-Uhlenbeck
#'
#' Reparameterise a phylogenetic OU process as a Gauss-Markov process over a
#' phylogeny.
#'
#' @inheritParams simulate_phylogenetic_ou
#'
#' @return A \eqn{2n - 1 \times 2} matrix of positive real values, where \eqn{n}
#'   the the number of terminal nodes on a bifurcating tree. The first column is
#'   the weight of the parent node to the expected value, while the second is
#'   the conditional standard deviation given the parent.
#' @export
#'
#'
#'
#' @examples
#' S <- 100
#' tr <- ape::rtree(S)
#' par <- reparameterise_phylogenetic_ou(tr)
reparameterise_phylogenetic_ou <- function(
  phy,
  heritable_amplitude = 1, length_scale = 2,
  environmental_amplitude = 0,
  perform_checks = TRUE
){
  if (perform_checks) {
    checkmate::assert_number(heritable_amplitude, lower = 0, finite = TRUE)
    checkmate::assert_number(length_scale, lower = 0, finite = TRUE)
    checkmate::assert_number(environmental_amplitude, lower = 0, finite = TRUE)
    if (is.null(phy$edge.length))
      stop("tree has no branch length")
    if (any(phy$edge.length < 0))
      stop("at least one branch length negative")
  }
  n <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- n + 1L
  par <- matrix(nrow = n + phy$Nnode, ncol = 2)
  par[ROOT, ] <- c(1, heritable_amplitude^2)
  des <- phy$edge[, 2]
  d <- phy$edge.length
  for (i in 1:N) {
    par[des[i], ] <- c(
      exp(- d[i] / length_scale),
      (heritable_amplitude^2)* (1 - exp(- 2 * d[i] / length_scale))  +
        (des[i] <= n) * environmental_amplitude^2
    )
  }
  par[, 2] <- sqrt(par[, 2])
  colnames(par) <- c("weight", "sd")
  par
}
