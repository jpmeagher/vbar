test_that("ordinal link function", {
  N <- 100
  gamma <- c(-Inf, 0, 0.5, 1, Inf)
  x <- rnorm(N)

  y <- sapply(x, function(z) sum(z > gamma))
  expect_true(
    all(x > gamma[y] & x < gamma[y+1])
  )

  y <- factor(y, ordered = TRUE, levels = 1:4)
  expect_equal(
    y,
    ordinal_link(x, gamma)
  )
})

test_that("ordinal inverse link", {
  N <- 100
  gamma <- c(-Inf, 0, 0.5, 1, Inf)
  x <- rnorm(N)
  y_num <- sapply(x, function(z) sum(z > gamma))
  expect_true(
    all(x > gamma[y_num] & x < gamma[y_num+1])
  )
  y <- factor(y_num, ordered = TRUE, levels = 1:4)

  x_map <- ordinal_inverse_link(
    y = y, cut_off_points = gamma,
    mu = rep(0, N)
  )
  expect_true(
    all(x_map > gamma[y_num] & x_map < gamma[y_num+1])
  )

  expect_true(
    all(sapply(1:4, function(i)all(x_map[y_num == i] == RcppTN::etn(.low = gamma[i], .high = gamma[i+1]))))
  )

  x_map <- ordinal_inverse_link(
    y = y, cut_off_points = gamma,
    mu = rep(0, N),
    return_expectation = FALSE
  )
  expect_true(
    all(x_map > gamma[y_num] & x_map < gamma[y_num+1])
  )

  x_map <- ordinal_inverse_link(
    y = y, cut_off_points = gamma,
    mu = rep(0, N),
    return_expectation = FALSE, random_seed = 1
  )
  expect_equal(
    x_map,
    ordinal_inverse_link(
      y = y, cut_off_points = gamma,
      mu = rep(0, N),
      return_expectation = FALSE, random_seed = 1
    )
  )
})

test_that("nominal link function", {
  N <- 100
  K <- nlevels(synthetic_traits$nom)
  X <- matrix(rnorm(N * K), nrow = N, ncol = K)
  categories <- factor(levels(synthetic_traits$nom))
  y <- apply(X, 1, which.max)

  max_check <- sapply(1:N, function(i)  all(X[i, y[i]] > X[i, -y[i]]))
  expect_true(
    all(max_check)
  )

  y <- categories[y]
  expect_equal(
    y,
    nominal_link(X, levels = categories)
  )
})

test_that("nominal trait normalising constant", {
  set.seed(101)
  n <- 10000
  K <- 10
  mu <- rnorm(K, sd = 0.5)
  x <- rnorm(K, mean = mu)
  y <- which.max(x)

  u <- rnorm(n)

  Z <- sapply(u + mu[y], function(x) x - mu[-y]) %>%
    pnorm() %>%
    apply(2, prod) %>%
    mean

  probs <- replicate(n, rnorm(K, mean = mu)) %>%
    apply(2, which.max) %>%
    table() %>%
    magrittr::divide_by(n)

  expect_equal(
    Z,
    unname(probs[y]),
    tolerance = 1 / (sqrt(n) * abs(unname(probs[y])) )
  )

  expect_equal(
    Z,
    nominal_probit_normalising_constant(
      i = y, mu = mu, n_samples = n
    ),
    tolerance = 1 / (sqrt(n) * Z)
  )
})

test_that("ordinal trait normalising constant", {
  gamma <- c(-Inf, 0, 0.5, 1, Inf)
  K <- length(gamma) - 1
  N <- 1000
  mu <- rep(0.5, N)
  X <- rnorm(N, mean = mu, sd = 0.5)
  y <- sapply(X, function(x) sum(x > gamma))
  expect_true(
    all(X > gamma[y] & X < gamma[y+1])
  )
  Z <- pnorm(gamma[y+1], mean = mu) - pnorm(gamma[y], mean = mu)

  y_fac <- factor(y, ordered = TRUE)
  expect_equal(
    ordinal_probit_normalising_constant(
      y = y_fac, mu = mu, cut_off_points = gamma,
      log_out = FALSE,
      perform_checks = TRUE
    ),
    Z
  )

  y_fac <- factor(y - 2, ordered = TRUE)
  expect_equal(
    exp(ordinal_probit_normalising_constant(
      y = y_fac, mu = mu, cut_off_points = gamma,
      perform_checks = TRUE
    )),
    Z
  )

  checkmate::expect_numeric(
    exp(ordinal_probit_normalising_constant(
      y = y_fac, mu = mu, cut_off_points = gamma,
      perform_checks = TRUE
    )),
    lower = 0, upper = 1, any.missing = FALSE
  )

})

test_that("auxiliary variable expectation",{
  set.seed(101)
  n <- 1e2
  K <- 10
  mu <- rnorm(K, sd = 0.5)
  x <- rnorm(K, mean = mu)
  i <- which.max(x)
  expected_auxiliary_1 <- expected_auxiliary_2 <- rep(NA, K)

  rs <- 101
  Z <- nominal_probit_normalising_constant(
    i = i, mu = mu, n_samples = n, random_seed = rs
  )

  set.seed(rs)
  u <- rnorm(n)
  z <- sapply(u, function(z){
    z + mu[i] - mu[-i]
  }) %>% t

  cdf <- pnorm(z)
  expected_auxiliary_1[i] <- mu[i] + mean(u * apply(cdf, 1, prod)) / Z

  den <- dnorm(z)
  tmp <- sapply(1:(K - 1), function(k){
    mean(den[, k] * apply(cdf[, -k], 1, prod))
  })
  expected_auxiliary_1[-i] <- expected_auxiliary_2[-i] <-  mu[-i] - (tmp / Z)
  expected_auxiliary_2[i] <- mu[i] + sum(mu[-i] - expected_auxiliary_2[-i])

  expect_equal(
    expected_auxiliary_2,
    nominal_probit_auxiliary_expectation(
      i = i, mu = mu, n_samples = n, random_seed = rs
    )
  )
})

test_that("nominal trait inverse link", {
  N <- 100
  K <- nlevels(synthetic_traits$nom)
  X <- matrix(rnorm(N * K), nrow = N, ncol = K)
  categories <- factor(levels(synthetic_traits$nom))
  y_num <- apply(X, 1, which.max)
  y <- categories[y_num]

  X_map <- nominal_inverse_link(
    y = y, mu = matrix(0, N, K),
    return_expectation = FALSE
  )

  max_check <- sapply(1:N, function(i)  all(X_map[i, y_num[i]] > X_map[i, -y_num[i]]))
  expect_true(
    all(max_check)
  )

  X_map <- nominal_inverse_link(
    y = y, mu = matrix(0, N, K),
    return_expectation = TRUE
  )
  max_check <- sapply(1:N, function(i)  all(X_map[i, y_num[i]] > X_map[i, -y_num[i]]))
  expect_true(
    all(max_check)
  )
})

test_that("initialise auxiliary traits", {
  P <- 4
  tn <- c("ord", "nom", "con", "fvt")
  tt <- factor(tn, levels = c("ord", "nom", "con", "fvt"))
  ind_mt <- list(
    ord = 2L, nom = 3L, con = 4L, fvt = 5:36
  )
  mt <- cbind(synthetic_traits[, 1:4], fvt = t(simplify2array(synthetic_traits$fvt)))
  N <- nrow(mt)
  K <- c(4, 3, NA, NA)
  ind_at <- list(
    ord = 1L, nom = 1 + 1:3, con = 5L, fvt = 5 + 1:32
  )
  g <- list(
    ord = function(x){
      ordinal_link(x, cut_off_points = gamma$ord)
    },
    nom = function(x) {
      nominal_link(x, levels = cat$nom)
    },
    con = function(x) x, fvt = function(x) exp(x)
  )
  g_inv <- list(
    ord = function(y){
      ordinal_inverse_link(
        y, cut_off_points = gamma$ord,
        mu = rep(0, N), return_expectation = FALSE
      )
    },
    nom = function(y){
      nominal_inverse_link(
        y, mu = matrix(0, N, K[2]),
        n_samples = 1000, return_expectation = FALSE
      )
    },
    con = function(y) data.matrix(y),
    fvt = function(y) log(data.matrix(y))
  )
  gamma <- list(
    ord = c(-Inf, 0, 1, 2, Inf),
    nom = NA, con = NA, fvt = NA
  )
  cat <- list(
    ord = NA,
    nom = factor(levels(mt$nom)),
    con = NA, fvt = NA
  )

  meta <- specify_manifest_trait_metadata(
    n_traits = P, trait_names = tn, trait_type = tt,
    trait_levels = K,
    manifest_trait_index = ind_mt, auxiliary_trait_index = ind_at,
    link_functions = g,
    inverse_link_functions = g_inv,
    cut_off_points = gamma, categories = cat,
    manifest_trait_df = mt,
    perform_checks = TRUE
  )

  X <- initialise_auxiliary_traits(
    n_traits = nrow(meta),
    manifest_trait_df = mt,
    trait_names = meta$trait_names,
    trait_type = meta$trait_type,
    trait_levels = meta$trait_levels,
    manifest_trait_index = meta$manifest_trait_index,
    auxiliary_trait_index = meta$auxiliary_trait_index,
    inverse_link_functions = meta$inverse_link_functions,
    auxiliary_traits = NULL,
    perform_checks = TRUE
  )

  for (i in 1:3) {
    expect_equal(
      g[[i]](X[, ind_at[[i]]]), mt[, ind_mt[[i]]]
    )
  }
  i <- 4
  expect_equal(
    g[[i]](X[, ind_at[[i]]]), data.matrix(mt[, ind_mt[[i]]]),
    ignore_attr = TRUE
  )

  expect_equal(
    initialise_auxiliary_traits(
      n_traits = nrow(meta),
      manifest_trait_df = mt,
      trait_names = meta$trait_names,
      trait_type = meta$trait_type,
      trait_levels = meta$trait_levels,
      manifest_trait_index = meta$manifest_trait_index,
      auxiliary_trait_index = meta$auxiliary_trait_index,
      inverse_link_functions = meta$inverse_link_functions,
      auxiliary_traits = X,
      perform_checks = TRUE
    ),
    X
  )
})

test_that("Auxiliary Precision Initialised", {
  P <- 4
  tn <- c("ord", "nom", "con", "fvt")
  tt <- factor(tn, levels = c("ord", "nom", "con", "fvt"))

  lambda <- c(1, 1, rgamma(2, shape = 1, rate = 0.01))
  names(lambda) <- tn

  prec <- initialise_precision(
    n_traits = P, trait_names = tn,  trait_type = tt,
    precision_prior_shape = 1, precision_prior_rate = 0.01,
    precision = NULL,
    perform_checks = TRUE
  )
  expect_equal(
    names(prec),
    tn
  )
  expect_equal(
    prec[tt %in% c("ord", "nom")],
    rep(1, sum(tt %in% c("ord", "nom"))),
    ignore_attr = TRUE
  )

  expect_equal(
    lambda,
    initialise_precision(
      n_traits = P, trait_names = tn,  trait_type = tt,
      precision_prior_shape = 1, precision_prior_rate = 0.01,
      precision = lambda,
      perform_checks = TRUE
    )
  )

  lambda["ord"] <- 1.1
  expect_error(
    initialise_precision(
      n_traits = P, trait_names = tn,  trait_type = tt,
      precision_prior_shape = 1, precision_prior_rate = 0.01,
      precision = lambda,
      perform_checks = TRUE
    )
  )
})

test_that("precision maps to vector", {
  P <- 4
  tn <- c("ord", "nom", "con", "fvt")
  tt <- factor(tn, levels = c("ord", "nom", "con", "fvt"))
  ind_at <- list(
    ord = 1L, nom = 1 + 1:3, con = 5L, fvt = 5 + 1:32
  )
  D_p <- sum(sapply(ind_at, length))

  lambda <- initialise_precision(
    n_traits = P, trait_names = tn,  trait_type = tt,
    precision_prior_shape = 1, precision_prior_rate = 0.01,
    precision = NULL,
    perform_checks = TRUE
  )

  lambda_vector <- map_precision_to_auxiliary_traits(
    precision = lambda, auxiliary_trait_index = ind_at
  )

  for(i in 1:P) {
    expect_equal(
      lambda[tn[i]],
      unique(lambda_vector[ind_at[[tn[i]]]]),
      ignore_attr = TRUE
    )
  }
})

test_that("Ordinal auxiliary trait elbo computed", {
  N <- 100
  L <- 4
  w <- rnorm(L)
  Z <- matrix(rnorm(L*N), N, L)
  w_outer <- (w %*% t(w)) + diag(runif(L))
  Z_outer <- lapply(
    1:N, function(i){
      (Z[i, ] %*% t(Z[i, ])) + diag(runif(L))
    }) %>%
    simplify2array()

  gamma <- c(-Inf, 0, sort(runif(2)), Inf)
  K <- length(gamma) -1
  X <- rnorm(N, mean = Z %*% w)
  y <- sapply(X, function(x) sum(x > gamma))
  expect_true(
    all(X > gamma[y] & X < gamma[y+1])
  )

  log_Z <- sum(log(pnorm(gamma[y + 1], mean = Z %*% w) - pnorm(gamma[y], mean = Z %*% w)))
  m_ex_2 <- sum((Z %*% w)^2)
  m_2_ex <- sum(sapply(
    1:N, function(i){
      sum(diag(w_outer %*% Z_outer[, , i]))
    }))

  y_fac <- factor(y - 2, ordered = TRUE, levels = -1:2)
  expect_equal(
    compute_ordinal_auxiliary_trait_elbo(
      y = y_fac, cut_off_points = gamma,
      loading_expectation = w, latent_trait_expectation = Z,
      loading_outer_expectation = w_outer, latent_trait_outer_expectation = Z_outer,
      perform_checks = TRUE
    ),
    log_Z + (0.5 *(m_ex_2 - m_2_ex))
  )
})

test_that("Nominal auxiliary trait elbo computed", {
  N <- 100
  L <- 4
  K <- 5
  W <- matrix(rnorm(L*K), K, L)
  Z <- matrix(rnorm(L*N), N, L)
  W_outer <- lapply(
    1:K, function(i){
      (W[i, ] %*% t(W[i, ])) + diag(runif(L))
    }) %>%
    simplify2array()
  Z_outer <- lapply(
    1:N, function(i){
      (Z[i, ] %*% t(Z[i, ])) + diag(runif(L))
    }) %>%
    simplify2array()
  M <- Z %*% t(W)
  X <- t(apply(M, 1, function(x) rnorm(K, mean = x)))
  y <- apply(X, 1, which.max)
  y_fac <- factor(y, ordered = FALSE)
  set.seed(101)
  log_Z <- sapply(
    1:N, function(j){
      nominal_probit_normalising_constant(
        y[j], mu = M[j, ], n_samples = 1000,
        random_seed = NULL,
        log_out = TRUE,
        perform_checks = TRUE
      )
    }
  )
  i <- 1
  m_2_ex <- sum(sapply(
    1:N, function(i){
      sum(sapply(
        1:K, function(j){
          sum(diag(W_outer[, , j] %*% Z_outer[, , i]))
        }
      ))
    }))

  expect_equal(
    compute_nominal_auxiliary_trait_elbo(
      y = y_fac, n_samples = 1000, random_seed = 101,
      loading_expectation = W, latent_trait_expectation = Z,
      loading_outer_expectation = W_outer, latent_trait_outer_expectation = Z_outer,
      perform_checks = TRUE
    ),
    sum(log_Z) + (0.5 * (sum(M^2) - m_2_ex)),
    tolerance = 0.001
  )
})

test_that("Continuous auxiliary trait elbo computed", {
  N <- 100
  L <- 4
  w <- rnorm(L)
  w_var <-  diag(runif(L))
  w_outer <- (w %*% t(w)) + w_var
  Z <- matrix(rnorm(L*N), N, L)
  Z_var <- diag(runif(L))
  Z_outer <- lapply(
    1:N, function(i){
      (Z[i, ] %*% t(Z[i, ])) + Z_var
    }) %>%
    simplify2array()
  M <- Z %*% w
  lambda <- 10
  X <- rnorm(N, mean = M, sd = sqrt(1 / lambda))

  M_2_ex <- sum(
    diag(
      apply(w_outer, c(1, 2), sum) %*% apply(Z_outer, c(1, 2), sum)
    )
  )

  expect_equal(
    compute_continuous_auxiliary_trait_elbo(
      auxiliary_trait = X,
      loading_expectation = w, latent_trait_expectation = Z, precision = lambda,
      loading_outer_expectation = w_outer, latent_trait_outer_expectation = Z_outer,
      perform_checks = TRUE
    ),
    - (0.5 * N * log( 2 * pi)) + (0.5 * N * log(lambda)) -
      (lambda * (sum(X^2) - (2 * c(t(X) %*% M)) + M_2_ex) / 2)
  )

  D_p <- 10
  W <- matrix(rnorm(L*D_p), D_p, L)
  W_var <-  lapply(1:D_p, function(i) diag(runif(L)))
  W_outer <- lapply(
    1:D_p, function(i){
      (W[i, ] %*% t(W[i, ])) + W_var[[i]]
    }) %>%
    simplify2array()

  M <- Z %*% t(W)
  X <-  M + matrix(rnorm(N*D_p, sd = sqrt(1 / lambda)), N, D_p)

  M_2_ex <- sum(
    diag(
      apply(w_outer, c(1, 2), sum) %*% apply(Z_outer, c(1, 2), sum)
    )
  )

  expect_equal(
    compute_continuous_auxiliary_trait_elbo(
      auxiliary_trait = X,
      loading_expectation = W, latent_trait_expectation = Z, precision = lambda,
      loading_outer_expectation = w_outer, latent_trait_outer_expectation = Z_outer,
      perform_checks = TRUE
    ),
    - (0.5 * N * D_p * log( 2 * pi)) + (0.5 * N * D_p * log(lambda)) -
      (lambda * (sum(X^2) - (2 * sum(diag(t(X) %*% M))) + M_2_ex) / 2)
  )
})

test_that("Auxiliary trait elbo computed",{
  P <- 4
  L <- 4
  tn <- c("ord", "nom", "con", "fvt")
  tt <- factor(tn, levels = c("ord", "nom", "con", "fvt"))
  ind_mt <- list(
    ord = 2L, nom = 3L, con = 4L, fvt = 5:36
  )
  mt <- cbind(synthetic_traits[, 1:4], fvt = t(simplify2array(synthetic_traits$fvt)))
  N <- nrow(mt)
  K <- c(4, 3, NA, NA)
  ind_at <- list(
    ord = 1L, nom = 1 + 1:3, con = 5L, fvt = 5 + 1:32
  )
  D_p <- sum(sapply(ind_at, length))

  cat <- list(
    ord = NA,
    nom = factor(levels(mt$nom)),
    con = NA, fvt = NA
  )
  gamma <- list(
    ord = c(-Inf, 0, 1, 2, Inf),
    nom = NA, con = NA, fvt = NA
  )

  g <- list(
    ord = ordinal_link, nom = nominal_link, con = function(x) x, fvt = function(x) exp(x)
  )
  g_inv <- list(
    ord = function(y){
      ordinal_inverse_link(
        y, cut_off_points = gamma$ord,
        mu = rep(0, N), return_expectation = FALSE
      )
    },
    nom = function(y){
      nominal_inverse_link(
        y, mu = matrix(0, N, length(cat$nom)),
        n_samples = 1000, return_expectation = FALSE
      )
    },
    con = function(y) y,
    fvt = function(y) log(data.matrix(y))
  )

  meta <- specify_manifest_trait_metadata(
    n_traits = P, trait_names = tn, trait_type = tt,
    trait_levels = K,
    manifest_trait_index = ind_mt, auxiliary_trait_index = ind_at,
    link_functions = g,
    inverse_link_functions = g_inv,
    cut_off_points = gamma, categories = cat,
    manifest_trait_df = mt,
    perform_checks = TRUE
  )

  C_w <- diag(D_p)
  x <- seq(0, 1, length.out = length(ind_at[[4]]))
  d <- abs(outer(x, x, "-"))
  ell <- 1 / (2 * pi)
  C_w[ind_at[[4]], ind_at[[4]]] <-  (exp_quad_kernel(d, 1, ell) + (1e-6 * diag(length(ind_at[[4]])))) / (1 + 1e-6)

  ph <- vbar::synthetic_trait_model_specification$phylogeny
  S <- length(ph$tip.label)

  plvm <- initialise_plvm(
    manifest_trait_df = mt, metadata = meta, phy = ph,
    L = L,
    loading_prior_correlation = C_w,
    auxiliary_traits = NULL,
    precision = NULL,
    ard_precision = NULL,
    ard_shape = 1, ard_rate = 1,
    loading = NULL, method = "random",
    within_taxon_amplitude = NULL,
    heritable_amplitude = NULL,
    length_scale = 2,
    perform_checks = TRUE
  )

  elbo <- compute_auxiliary_trait_elbo(
    manifest_trait_df = mt, metadata = meta,
    auxiliary_traits = plvm$auxiliary_traits,
    loading_expectation = plvm$loading_expectation,
    latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
    precision = plvm$precision,
    loading_outer_expectation = plvm$loading_row_outer_product_expectation,
    latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
    n_samples = 1000, random_seed = 101,
    perform_checks = TRUE
  )

  expect_equal(
    elbo,
    compute_ordinal_auxiliary_trait_elbo(
      y = mt$ord, cut_off_points = gamma[[1]],
      loading_expectation = plvm$loading_expectation[meta$auxiliary_trait_index[[1]], ],
      latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
      loading_outer_expectation = plvm$loading_row_outer_product_expectation[, , meta$auxiliary_trait_index[[1]]],
      latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation
    ) +
      compute_nominal_auxiliary_trait_elbo(
        y = mt$nom, n_samples = 1000, random_seed = 101,
        loading_expectation = plvm$loading_expectation[meta$auxiliary_trait_index[[2]], ],
        latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
        loading_outer_expectation = plvm$loading_row_outer_product_expectation[, , meta$auxiliary_trait_index[[2]]],
        latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation
      ) +
      compute_continuous_auxiliary_trait_elbo(
        auxiliary_trait = plvm$auxiliary_traits[, meta$auxiliary_trait_index[[3]]],
        loading_expectation = plvm$loading_expectation[meta$auxiliary_trait_index[[3]], ],
        latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
        precision = plvm$precision[3],
        loading_outer_expectation = plvm$loading_row_outer_product_expectation[, , meta$auxiliary_trait_index[[3]]],
        latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
        perform_checks = TRUE) +
      compute_continuous_auxiliary_trait_elbo(
        auxiliary_trait = plvm$auxiliary_traits[, meta$auxiliary_trait_index[[4]]],
        loading_expectation = plvm$loading_expectation[meta$auxiliary_trait_index[[4]], ],
        latent_trait_expectation = plvm$individual_specific_latent_trait_expectation,
        precision = plvm$precision[4],
        loading_outer_expectation = plvm$loading_row_outer_product_expectation[, , meta$auxiliary_trait_index[[4]]],
        latent_trait_outer_expectation = plvm$individual_specific_latent_trait_outer_product_expectation,
        perform_checks = TRUE)
  )
})

test_that("Auxiliary Traits Updated", {
  P <- 4
  L <- 4
  tn <- c("ord", "nom", "con", "fvt")
  tt <- factor(tn, levels = c("ord", "nom", "con", "fvt"))
  ind_mt <- list(
    ord = 2L, nom = 3L, con = 4L, fvt = 5:36
  )
  mt <- cbind(synthetic_traits[, 1:4], fvt = t(simplify2array(synthetic_traits$fvt)))
  N <- nrow(mt)
  K <- c(4, 3, NA, NA)
  ind_at <- list(
    ord = 1L, nom = 1 + 1:3, con = 5L, fvt = 5 + 1:32
  )
  D_p <- sum(sapply(ind_at, length))

  cat <- list(
    ord = NA,
    nom = factor(levels(mt$nom)),
    con = NA, fvt = NA
  )
  gamma <- list(
    ord = c(-Inf, 0, 1, 2, Inf),
    nom = NA, con = NA, fvt = NA
  )

  g <- list(
    ord = ordinal_link, nom = nominal_link, con = function(x) x, fvt = function(x) exp(x)
  )
  g_inv <- list(
    ord = function(y){
      ordinal_inverse_link(
        y, cut_off_points = gamma$ord,
        mu = rep(0, N), return_expectation = FALSE
      )
    },
    nom = function(y){
      nominal_inverse_link(
        y, mu = matrix(0, N, length(cat$nom)),
        n_samples = 1000, return_expectation = FALSE
      )
    },
    con = function(y) y,
    fvt = function(y) log(data.matrix(y))
  )

  meta <- specify_manifest_trait_metadata(
    n_traits = P, trait_names = tn, trait_type = tt,
    trait_levels = K,
    manifest_trait_index = ind_mt, auxiliary_trait_index = ind_at,
    link_functions = g,
    inverse_link_functions = g_inv,
    cut_off_points = gamma, categories = cat,
    manifest_trait_df = mt,
    perform_checks = TRUE
  )

  C_w <- diag(D_p)
  x <- seq(0, 1, length.out = length(ind_at[[4]]))
  d <- abs(outer(x, x, "-"))
  ell <- 1 / (2 * pi)
  C_w[ind_at[[4]], ind_at[[4]]] <-  (exp_quad_kernel(d, 1, ell) + (1e-6 * diag(length(ind_at[[4]])))) / (1 + 1e-6)

  ph <- vbar::synthetic_trait_model_specification$phylogeny
  S <- length(ph$tip.label)

  plvm <- initialise_plvm(
    manifest_trait_df = mt, metadata = meta, phy = ph,
    L = L,
    loading_prior_correlation = C_w,
    auxiliary_traits = NULL,
    precision = NULL,
    ard_precision = NULL,
    ard_shape = 1, ard_rate = 1,
    loading = NULL, method = "random",
    within_taxon_amplitude = NULL,
    heritable_amplitude = NULL,
    length_scale = 2,
    perform_checks = TRUE
  )

  X <- update_discrete_auxiliary_traits(
    manifest_trait_df = mt, metadata = meta,
    auxiliary_traits = plvm$auxiliary_traits,
    loading_expectation = plvm$loading_expectation,
    latent_trait_expectation = plvm$individual_specific_latent_trait_expectation
  )

  expect_equal(
    X[, -(1:4)], plvm$auxiliary_traits[, -(1:4)]
  )

  expect_equal(
    sapply(X[, 1], function(x) sum(x > meta$cut_off_points[[1]])),
    as.numeric(mt$ord)
  )
  expect_equal(
    sapply(plvm$auxiliary_traits[, 1], function(x) sum(x > meta$cut_off_points[[1]])),
    as.numeric(mt$ord)
  )

  expect_equal(
    apply(X[, 2:4], 1, which.max),
    as.numeric(mt$nom)
  )
  expect_equal(
    apply(plvm$auxiliary_traits[, 2:4], 1, which.max),
    as.numeric(mt$nom)
  )
})
