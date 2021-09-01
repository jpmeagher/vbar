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

test_that("ordianl inverse link", {
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

test_that("auxiliary variavle expectation",{
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
