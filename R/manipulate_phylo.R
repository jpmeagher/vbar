#' Rescale a Phylogenetic tree
#'
#' Rescale branch lengths for a \eqn{phylo} object to fix the maximum distance
#' from the root to tip of a phylogenetic tree.
#'
#' @inheritParams simulate_phylogenetic_ou
#' @param max_dist A positive, real scalar. The maximum distance along the
#'   rescaled phylogeny from root to tip.
#'
#' @return An object of class "\eqn{phylo}" where the maximum distance from the
#'   root to tip is \eqn{max_dist}
#' @export
#'
#' @examples
#' tr <- ape::rtree(100)
#' tr1 <- scale_phylo(tr)
scale_phylo <- function(phy, max_dist = 1){
  checkmate::expect_number(max_dist, lower = 0)
  m <- max(phytools::nodeHeights(phy))
  phy$edge.length <- max_dist * (phy$edge.length / m)
  phy
}
