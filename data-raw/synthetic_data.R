S <- 2^7
set.seed(98)
phy <- ape::rcoal(S) %>%
  scale_phylo(max_dist = 1)

## Model Specification

P <- 4
N_s <- rep(3, S)
N <- sum(N_s)
K_ord <- 4
K_nom <- 3

## Set Loading

auxiliary_trait_precision <- c(ord = 1, nom = 1, con = 10, fvt = 10000)
L <- 3
# FVT Loading
D_fvt <- 2^5
x_fvt <- seq(0, 1, length.out = D_fvt)
loading_fvt_means <- c(.2, .5, .75)
loading_fvt_sds <- c(.06, .12, .07)
loading_fvt <- sapply(
  seq_along(loading_fvt_means),
  function(i) dnorm(x_fvt, mean = loading_fvt_means[i], sd = loading_fvt_sds[i])
)
loading_fvt <- loading_fvt %*% diag(1 / c(5, 4, 6))
# Continuous Trait Loading
loading_con <- c(0.5, 0.25, -0.25)
# Nominal Trait loading
loading_nom <- matrix(c(1, 0, 0, 0, -0.5, 0, 0, 0, 1), byrow = TRUE, nrow = K_nom, ncol = L)
# Ordinal Trait Loading
loading_ord <- c(0, 1, 0.5)
# Full Loading
loading <- rbind(loading_ord, loading_nom, loading_con, loading_fvt)
rownames(loading) <- c(
  "ord",
  paste("nom", 1:K_nom, sep = "."),
  "con",
  paste("fvt", 1:D_fvt, sep = ".")
)
colnames(loading) <- paste("loading", 1:L, sep = ".")

## Set Latent Traits

tau <- sqrt(c(0.001, 0.05, 0.1))
h <- sqrt(c(0.95, 0.66, 0.25))

set.seed(99)
taxon_specific_latent_traits <- sapply(
  h, function(x) {
    simulate_phylogenetic_ou(
      phy = phy,
      heritable_amplitude = x, length_scale = 2,
      environmental_amplitude = sqrt(1 - x^2),
      internal = TRUE
    )
  }
)

taxon_id <- sapply(1:S, function(i) rep(phy$tip.label[i], N_s[i])) %>%
  c()

individual_specific_latent_traits <- sapply(
  1:L, function(l){
    sqrt(1 - tau[l]^2) * taxon_specific_latent_traits[taxon_id, l] +
      tau[l] * rnorm(N)
  })

## Set Auxiliary Traits

D_prime <- nrow(loading)
auxiliary_traits_expectation <- individual_specific_latent_traits %*% t(loading)
auxiliary_trait_precision_vec <- c(
  auxiliary_trait_precision[1],
  rep(auxiliary_trait_precision[2], K_nom),
  auxiliary_trait_precision[3],
  rep(auxiliary_trait_precision[4], D_fvt)
)

set.seed(100)
auxiliary_traits <- auxiliary_traits_expectation +
  t(
    sapply(1:N, function(i){
      rnorm(D_prime, sd = sqrt(1 / auxiliary_trait_precision_vec))
    }
    )
  )

## Set Manifest Traits

cut_off_points <- c(-Inf, 0, 0.5, 1, Inf)

auxiliary_index_ord <- 1
auxiliary_index_nom <- max(auxiliary_index_ord) + (1:K_nom)
auxiliary_index_con <- max(auxiliary_index_nom) + 1
auxiliary_index_fvt <- max(auxiliary_index_con) + 1:D_fvt

synthetic_traits <- data.frame(
  taxon_id = taxon_id,
  ord = sapply(auxiliary_traits[, auxiliary_index_ord], function(x){
    x > cut_off_points
  }) %>%
    colSums() %>%
    factor(ordered = TRUE),
  nom = apply(auxiliary_traits[, auxiliary_index_nom], 1, which.max) %>%
    factor(ordered = FALSE),
  con = auxiliary_traits[, auxiliary_index_con]
)
synthetic_traits$fvt <- lapply(1:N, function(n) exp(unname(auxiliary_traits[n, auxiliary_index_fvt])))
usethis::use_data(synthetic_traits, overwrite = TRUE)

synthetic_trait_model_specification <- list(
  phylogeny = phy,
  auxiliary_traits = auxiliary_traits,
  loading = loading,
  auxiliary_trait_precision = auxiliary_trait_precision,
  auxiliary_trait_precision_vec = auxiliary_trait_precision_vec,
  individual_specific_latent_traits = individual_specific_latent_traits,
  taxon_specific_latent_traits = taxon_specific_latent_traits,
  within_taxon_sd = tau,
  heritability = h,
  ordinal_trait_cut_off = cut_off_points
)

usethis::use_data(synthetic_trait_model_specification, overwrite = TRUE)
