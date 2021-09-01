test_that("Ornstein-Uhlenbeck kernel", {
  x <- 1:10
  d <- abs(outer(x, x, "-"))
  alpha <-  5
  ell <- 2
  K <- ou_kernel(d, amplitude = alpha, length_scale = ell)

  expect_equal(
    K, alpha^2 * exp(- d / ell)
  )

  expect_equal(
    dim(d), dim(K)
  )
  expect_true(
    all(K > 0)
  )
  expect_true(
    all(K <= alpha^2)
  )

  expect_error(
    ou_kernel(-1, amplitude = 1, length_scale = 2)
  )
  expect_error(
    ou_kernel(-Inf, amplitude = 1, length_scale = 2)
  )
  expect_error(
    ou_kernel(0, amplitude = -1, length_scale = 2)
  )
  expect_error(
    ou_kernel(0, amplitude = Inf, length_scale = 2)
  )
  expect_error(
    ou_kernel(0, amplitude = 1, length_scale = -2)
  )
  expect_error(
    ou_kernel(0, amplitude = 1, length_scale = Inf)
  )
})

test_that("Exponentiated Quadratic kernel", {
  x <- 1:10
  d <- abs(outer(x, x, "-"))
  alpha <-  5
  ell <- 2
  K <- exp_quad_kernel(d, amplitude = alpha, length_scale = ell)

  expect_equal(
    K, (alpha^2 )* exp(- (d^2 )/ (2 * (ell^2)))
  )

  expect_equal(
    dim(d), dim(K)
  )
  expect_true(
    all(K > 0)
  )
  expect_true(
    all(K <= alpha^2)
  )

  expect_error(
    exp_quad_kernel(-1, amplitude = 1, length_scale = 2)
  )
  expect_error(
    exp_quad_kernel(-Inf, amplitude = 1, length_scale = 2)
  )
  expect_error(
    exp_quad_kernel(0, amplitude = -1, length_scale = 2)
  )
  expect_error(
    exp_quad_kernel(0, amplitude = Inf, length_scale = 2)
  )
  expect_error(
    exp_quad_kernel(0, amplitude = 1, length_scale = -2)
  )
  expect_error(
    exp_quad_kernel(0, amplitude = 1, length_scale = Inf)
  )
})

