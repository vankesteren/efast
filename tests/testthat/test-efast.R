context("EFAST testing")

test_that("EFAST fitting works", {
  # get test data
  test_data <- lavaan::HolzingerSwineford1939[,7:15]
  res_struct <- list(
    c("x4", "x7"),
    c("x5", "x9")
  )
  test_efast <- efast(test_data, 3, res_struct)

  expect_true(is_efast(test_efast))
  expect_true(inherits(efast_loadings(test_efast), "loadings"))
})

test_that("EFAST fitting with a covariance matrix works", {
  # get test data
  test_S <- cov(lavaan::HolzingerSwineford1939[,7:15])
  test_N <- nrow(lavaan::HolzingerSwineford1939)
  res_struct <- list(
    c("x4", "x7"),
    c("x5", "x9")
  )
  test_efast <- efast(test_S, 3, res_struct, sample.nobs = test_N)

  expect_true(is_efast(test_efast))
  expect_true(inherits(efast_loadings(test_efast), "loadings"))
})

test_that("EFAST fitting with a single latent variable works", {
  # get test data
  test_S <- cov(lavaan::HolzingerSwineford1939[,7:15])
  test_N <- nrow(lavaan::HolzingerSwineford1939)
  res_struct <- list(
    c("x4", "x7"),
    c("x5", "x9")
  )
  test_efast <- efast(test_S, 1, res_struct, sample.nobs = test_N)

  expect_true(is_efast(test_efast))
  expect_true(inherits(efast_loadings(test_efast), "loadings"))
})

test_that("EFA fitting works", {
  fit_efa <- efast_efa(data = lavaan::HolzingerSwineford1939[,7:15], M = 3)
  expect_true(is_efast_efa(fit_efa))
  expect_true(inherits(efast_loadings(fit_efa), "loadings"))
})

test_that("EFAST fitting on big data works", {
  set.seed(45)
  test_dat <- matrix(rnorm(5000), 100)
  test_dat <- as.data.frame(
    test_dat %*% chol(toeplitz(c(1, rep(0, 24), 0.5, rep(0, 24))))
  )
  fit_efast_big <- efast_hemi(
    data = test_dat,
    M = 3,
    lh_idx = 1:25,
    rh_idx = 26:50,
    constrain = FALSE
  )
  expect_true(is_efast_hemi(fit_efast_big))
})

test_that("Data generation works and efast thereon works", {
  set.seed(45)
  simdat <- simulate_efast()
  fit_sim <- efast_hemi(simdat, M = 4, 1:17, 18:34)
  expect_true(is_efast_hemi(fit_sim))
})

test_that("Methods on efast objects work", {
  set.seed(45)
  simdat <- simulate_efast()
  fit_sim <- efast_hemi(simdat, M = 4, 1:17, 18:34)
  expect_true(is_efast(fit_sim))
  li_tab <- lateralization(fit_sim)
  expect_equal(
    round(li_tab[,2], 3),
    c(0.314, 0.299, 0.304, 0.321, 0.316, 0.282, 0.321, 0.274, 0.218,
      0.219, 0.199, 0.217, 0.218, 0.200, 0.207, 0.221, 0.218)
  )
  comps <- decomposition(fit_sim)
  expct <- (comps$factor + comps$residual + comps$structure)
  expect_true(all(abs(comps$observed - expct) < 0.3))
})
