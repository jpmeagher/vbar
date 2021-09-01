#' Check trait type specification
#'
#' Check that specified trait type metadata matches the dataframe of manifest
#' traits. This provides a non-exhaustive series of checks.
#'
#' @param manifest_trait Variable(s) from a daframe. The manifest trait
#'   recordings to be checked.
#' @param name A character scalar. The name of the trait to be checked.
#' @param type A caracter scalar. One of \eqn{c("ord", "nom", "con", "fvt")}.
#'   The type of trait to be checked.
#' @param levels A scalar. The number of levels associated with a discrete manifest
#'   trait. NA when traits are continuous.
#'
#' @return If the conditions are not met then returns an error message.
check_trait_type <- function(
  manifest_trait, name, type, levels
){
  checkmate::assert_choice(
    as.character(type), c("ord", "nom", "con", "fvt"),
  )
  if (type == "ord") {
    checkmate::assert_factor(
      manifest_trait, ordered = TRUE, any.missing = FALSE, n.levels = levels,
      .var.name = name
    )
  }
  if (type == "nom") {
    checkmate::assert_factor(
      manifest_trait, ordered = FALSE, any.missing = FALSE, empty.levels.ok = FALSE, n.levels = levels,
      .var.name = name
    )
  }
  if (type == "con") {
    checkmate::assert_numeric(
      manifest_trait, any.missing = FALSE, .var.name = name
    )
    if (checkmate::test_integerish(manifest_trait)) {
      stop(paste("Error. assertion on '", name, "' failed. Must not be of type 'integerish'"))
    }
    checkmate::assert_scalar_na(levels, .var.name = paste("number of levels associated with", type))
  }
  if (type == "fvt") {
    checkmate::assert_numeric(
      unlist(manifest_trait), any.missing = FALSE, .var.name = name
    )
    if (checkmate::test_integerish(manifest_trait)) {
      stop(paste("Error. assertion on '", name, "' failed. Must not be of type 'integerish'"))
    }
    if (checkmate::test_vector(manifest_trait, strict = TRUE)) {
      stop(paste("Error. assertion on '", name, "' failed. Must not be of type 'vector'"))
    }
    checkmate::assert_scalar_na(levels, .var.name= paste("number of levels associated with", type))
  }
}

#' Compare Metadata to Manifest Traits
#'
#' Compare the supplied meta data to that implicit in a data frame of manifest
#' traits.
#'
#' @param n_traits A natural number. The number of manifest traits.
#' @param manifest_trait_df A data frame. Observed manifest traits.
#' @param trait_names A P-dimensional character vector. The manifest trait
#'   names.
#' @param trait_type A P-dimensional vector of unordered factors belonging to
#'   one of 4 levels. The trait type corresponding to each of the named manifest
#'   traits. Ordinal traits are denoted \eqn{ord}, nominal are \eqn{nom},
#'   scalar-valued continuous are \eqn{con}, and function-valued traits are
#'   \eqn{fvt}.
#' @param trait_levels A P-dimensional vector. The number of levels associated with each
#'   manifest trait. NA for continuous traits.
#' @param manifest_trait_index A P-dimensional list. columns of
#'   \eqn{manifest_trait_df} associated with each trait. List elements will be
#'   scalar valued for ordinal, nominal, and scalar-valued continuous traits.
#'   Function-valued traits are multi-variate and indexed by a vector of values.
#'
#' @return Reports error messages if metadata and dataframe do not match.
compare_metadata_to_df <- function(
  n_traits,
  manifest_trait_df,
  trait_names,
  trait_type,
  trait_levels,
  manifest_trait_index
){
  if (is.null(manifest_trait_df)) return(invisible(NULL))

  do_manifest_traits_match_trait_names_and_indices <-
    all(
      sapply(
        1:n_traits,
        function(i){
          tmp <- which(grepl(trait_names[i], colnames(manifest_trait_df))) == manifest_trait_index[[i]]
          all(tmp) & (length(tmp) > 0)
        }
      )
    )
  checkmate::assert_true(
    do_manifest_traits_match_trait_names_and_indices
  )
  lapply(
    1:n_traits, function(i){
      check_trait_type(
        manifest_trait = manifest_trait_df[, manifest_trait_index[[i]]],
        name = trait_names[i], type = trait_type[i], levels = trait_levels[i]
      )
    }
  )
  if (sum(sapply(manifest_trait_index, length)) < (ncol(manifest_trait_df) - 1)) {
    warning("Some traits in 'manifest_trait_df' may not be covered by 'manifest_trait_index'.")
  }

  invisible(NULL)
}

#' Specify Manifest Trait Metadata
#'
#' Converts metadata relating to a data frame of manifest traits into a
#' dataframe. Ensures that the meta data is correctly formatted, allowing CAVI
#' for comparative analysis and ancestral reconstruction,
#'
#' @inheritParams compare_metadata_to_df
#' @param auxiliary_trait_index A P-dimensional list. Columns of the auxiliary
#'   trait matrix to be associated with each trait. List elements will be scalar
#'   valued for ordinal, nominal, and scalar-valued continuous traits.
#'   Function-valued traits are multi-variate and indexed by a vector of values.
#' @param link_functions A P-dimensional list. Functions mapping auxiliary
#'   traits to manifest traits. For traits taking continuous values these must
#'   be specified by the user and can be the identity function.
#' @param inverse_link_functions A P-dimensional list. Functions mapping
#'   manifest traits to auxiliary traits. For traits taking continuous values
#'   these must be specified by the user and can be the identity function.
#' @param cut_off_points A P-dimensional list. The set of cut-off points for
#'   ordinal traits. \eqn{NA} for all other traits.
#' @param categories A P-dimensional list. The set of categories for nominal
#'   traits, stored as an unordered factor. \eqn{NA} for all other traits.
#' @param manifest_trait_df The data frame of manifest traits. If included
#'   allows for a non-exhaustive series of checks comparing metadata provided to
#'   that implicit in the manifest traits.
#' @inheritParams  ou_kernel
#'
#' @return A data frame containing required metadata for each manifest trait.
#'
#' @export
#'
#' @examples
#' P <- 4
#' tn <- c("ord", "nom", "con", "fvt")
#' tt <- factor(tn, levels = c("ord", "nom", "con", "fvt"))
#' ind_mt <- list(
#'   ord = 2L, nom = 3L, con = 4L, fvt = 5:36
#' )
#' mt <- cbind(synthetic_traits[, 1:4], fvt = t(simplify2array(synthetic_traits$fvt)))
#' K <- c(4, 3, NA, NA)
#' ind_at <- list(
#'   ord = 1L, nom = 1 + 1:3, con = 5L, fvt = 5 + 1:32
#' )
#' g <- list(
#'   ord = ordinal_link,
#'   nom = nominal_link,
#'   con = function(x) x,
#'   fvt = function(x) exp(x)
#' )
#' g_inv <- list(
#'   ord = ordinal_inverse_link,
#'   nom = nominal_inverse_link,
#'   con = function(y) y,
#'   fvt = function(y) log(y)
#' )
#' gamma <- list(
#'   ord = c(-Inf, 0, 1, 2, Inf),
#'   nom = NA, con = NA, fvt = NA
#' )
#' cat <- list(
#'   ord = NA,
#'   nom = factor(levels(mt$nom)),
#'   con = NA, fvt = NA
#' )
#'
#' meta <- specify_manifest_trait_metadata(
#'   n_traits = P, trait_names = tn, trait_type = tt,
#'   trait_levels = K,
#'   manifest_trait_index = ind_mt, auxiliary_trait_index = ind_at,
#'   link_functions = g,
#'   inverse_link_functions = g_inv,
#'   cut_off_points = gamma, categories = cat,
#'   manifest_trait_df = mt,
#'   perform_checks = TRUE
#' )
#' checkmate::expect_data_frame(meta)
specify_manifest_trait_metadata <- function(
  n_traits,
  trait_names, trait_type, trait_levels,
  manifest_trait_index,
  auxiliary_trait_index,
  link_functions,
  inverse_link_functions,
  cut_off_points, categories,
  manifest_trait_df = NULL,
  perform_checks = TRUE
){
  if (perform_checks) {
    checkmate::assert_character(
      trait_names, any.missing = FALSE, len = n_traits, unique = TRUE
    )
    checkmate::assert_factor(
      trait_type, levels = c("ord", "nom", "con", "fvt"),  ordered = FALSE,
      empty.levels.ok = TRUE, any.missing = FALSE, len = n_traits
    )
    checkmate::assert_list(
      manifest_trait_index, type = "numeric", any.missing = FALSE,
      len = n_traits, names = "named"
    )
    checkmate::assert_list(
      auxiliary_trait_index, type = "numeric", any.missing = FALSE,
      len = n_traits, names = "named"
    )
    checkmate::assert_true(
      all(names(manifest_trait_index) == trait_names)
    )
    checkmate::assert_set_equal(names(manifest_trait_index), names(auxiliary_trait_index))
    checkmate::assert_data_frame(
      manifest_trait_df, any.missing = FALSE,
      min.cols = sum(unlist(lapply(manifest_trait_index, length))),
      col.names = "named"
    )
    checkmate::assert_list(
      link_functions, types = "function", any.missing = FALSE, len = n_traits
    )
    checkmate::assert_list(
      inverse_link_functions, types = "function", any.missing = FALSE, len = n_traits
    )
    checkmate::assert_list(
      cut_off_points, len = n_traits
    )
    checkmate::assert_list(
      categories, len = n_traits
    )
    compare_metadata_to_df(
      n_traits = n_traits, manifest_trait_df = manifest_trait_df,
      trait_names = trait_names, trait_type = trait_type,
      trait_levels = trait_levels,
      manifest_trait_index = manifest_trait_index
    )
    manifest_D_prime <- sum(
      sapply(1:n_traits, function(i){
        if (trait_type[i] == "nom") {
          nlevels(manifest_trait_df[, manifest_trait_index[[i]]])
        } else {
          length(manifest_trait_index[[i]])
        }
      })
    )
    auxiliary_D_prime <- length(unique(unlist(auxiliary_trait_index)))
    checkmate::assert_true(
      manifest_D_prime == auxiliary_D_prime,
      .var.name = "Auxiliary trait index matches to manifest trait dimensions"
    )
    checkmate::assert_true(
      all(unlist(auxiliary_trait_index) %in% 1:auxiliary_D_prime),
      .var.name = "Auxiliary trait index runs from 1 to D'"
    )
  }
  out <- data.frame(
    trait_names = trait_names,
    trait_type = trait_type,
    trait_levels = trait_levels
  )
  out$manifest_trait_index <- manifest_trait_index
  out$auxiliary_trait_index <- auxiliary_trait_index
  out$link_functions <- link_functions
  out$inverse_link_functions <- inverse_link_functions
  out$cut_off_points <- cut_off_points
  out$categories <- categories

  out
}
