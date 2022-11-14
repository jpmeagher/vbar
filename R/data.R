#' Synthetic Phenotypic Traits
#'
#' A set of synthetic manifest traits generated sampled from a Phylogenetic
#' Latent Variable Model (PLVM). Used to illustrate the approach to ancestral
#' reconstruction offered by Variational Inference for the PLVM.
#'
#' @format A data frame with 384 observations of 36 variables:
#' \describe{
#'   \item{taxon_id}{The taxon to which the individual belongs.}
#'   \item{ord}{An
#'   ordinal trait for which K = 4}
#'   \item{nom}{A nominal trait for which K = 3}
#'   \item{con}{A scalar-valued continuous trait.}
#'   \item{fvt}{A function-valued trait lying on the positive real
#'   numbers, observed at 32 points registered along the continuum.}
#'}
"synthetic_traits"

#' Model Specification for Synthetic Traits
#'
#' The Phylogenetic Latent Variable Model (PLVM) specified in the package
#' README. Used to illustrate the approach to ancestral reconstruction offered
#' by Variational Inference for the PLVM.
#'
#' @format A list of length 9:
#' \describe{
#'   \item{phylogeny}{An object of type \eqn{phylo} with \eqn{S = 128} terminal
#'    taxa and \eqn{S-1} internal taxa}
#'    \item{auxiliary_traits}{A named 384x37 matrix of real numbers. Auxiliary
#'    traits mapping to synthetic manifest traits.}
#'    \item{loading}{A named 37x4 matrix of real numbers}
#'    \item{auxiliary_trait_precision}{A 4-dimensional vector of positive real
#'    numbers. Precision with which auxiliary traits are observed.}
#'    \item{auxiliary_trait_precision_vec}{A D'-dimensional vector of positive real
#'    numbers. Precision for each auxiliary trait.}
#'    \item{individual_specific_latent_traits}{A 384x3 matrix of real numbers.
#'    The latent trait associated with each observed manifest trait.}
#'    \item{taxon_specific_latent_traits}{A 255x3 matrix of real numbers.
#'    The latent trait associated with each taxa.}
#'    \item{within_taxon_sd}{A 3 dimensional vector of real values on the
#'    interval (0, 1). Standard deviation of the within-taxon variation.}
#'    \item{heritability}{A 3 dimensional vector of real values on the
#'    interval (0, 1). The heritability of latent traits over the phylogeny.}
#'    \item{ordinal_trait_cut_off}{A 5-dimensional vector of increasing values.
#'    Cut-off points for the ordinal trait mapping}
#'}
"synthetic_trait_model_specification"

#' Bat Phylogeny
#'
#' A phylogenetic tree for 22 species of Mexican bat transcribed from the
#' maximum a posteriori phylogeny inferred by Amador et al. (2018) (see link
#' below).
#'
#' @source{Bat Systematics in the Light of Unconstrained Analyses of a
#' Comprehensive Molecular Supermatrix -
#' https://link.springer.com/article/10.1007/s10914-016-9363-8}
"bat_phylogeny"

#' Bat scholocation call recordings
#'
#' A set of 1805 echolocation call recordings collected from 449 individual bats
#' across 22 species and 5 families. Calls were recorded at a sampling rate of
#' 500 kHz. The original data set is described by Stathopolous et al. (2017)
#' (see link below). Additional features have been extracted from echolocation
#' calls allow comparative analysis. These feature are based on a harmonic model
#' for bat echolocation calls and the feature extraction procedure is described
#' in Meagher (2020) (see second link below)
#'
#' @source{Bat echolocation call identification for biodiversity monitoring: a probabilistic approach - https://doi.org/10.1111/rssc.12217}
#' @source{Bayesian Ancestral Reconstruction for Bat Echolocation - https://ethos.bl.uk/OrderDetails.do?uin=uk.bl.ethos.828108}
#'
#' @format A data frame with 1816 observations of 4 variables:
#' \describe{
#' \item{bat}{The individual bat identifier.}
#' \item{species}{The bat species identifier.}
#' \item{family}{The species family identifier.}
#' \item{call_recording}{A list of calls recorded at a sampling rate of 500 kHz.}
#' \item{fundamental_frequency_contour}{A list of estimates for the fundamental
#' frequency taken from the bat echolocation call recording and the associated time
#' for the estimate within the call recording.}
#' \item{duration}{The time duration of the echolocation call. Reported in
#' milliseconds.}
#' \item{harmonic_order}{The number of harmonics in the echolocation call.}
#' \item{dominant_harmonic}{The harmonic carrying the most energy.}
#' }
"bat_echolocation_calls"
