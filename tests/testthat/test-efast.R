context("EFAST testing")

test_that("EFAST fitting works", {
  fit_efast <- efast_hemi(
    data      = lavaan::HolzingerSwineford1939,
    M         = 3,
    lh_idx    = 7:10,
    rh_idx    = 11:14,
    roi_names = c("a", "b", "c", "d"),
    constrain = TRUE
  )
  expect_true(is_efast(fit_efast))
  expect_true(inherits(efast_loadings(fit_efast), "loadings"))
})



test_that("EFA fitting works", {
  fit_efa <- efa_esem(data = lavaan::HolzingerSwineford1939[,7:14], M = 3)
  expect_true(is_efa(fit_efa))
  expect_true(inherits(efast_loadings(fit_efa), "loadings"))
})


test_that("EFAST fitting on big data works and methods work", {
  set.seed(45)
  test_dat <- as.data.frame(matrix(rnorm(5000), 100))
  fit_efast_big <- efast_hemi(
    data = test_dat,
    M = 3,
    lh_idx = 1:25,
    rh_idx = 26:50,
    constrain = FALSE
  )
  expect_true(is_efast(fit_efast_big))
  li_tab <- lateralization(fit_efast_big)
  expect_equal(
    round(li_tab[,2], 3),
    c(0.967, 0.997, 1, 1, 1, 0.44, 1, 0.839, 1, 1, 1, 0.943, 1, 0.888, 0.891,
      1, 1, 0.984, 0.837, 0.749, 0.937, 1, 1, 1, 0.95)
  )
  comps <- decomposition(fit_efast_big)
  expct <- (comps$factor + comps$residual + comps$structure)
  expect_true(all(abs(comps$observed - expct) < 0.3))
})

test_that("Data generation works and efast thereon works", {
  simdat <- simulate_efast()
  fit_sim <- efast_hemi(simdat, M = 4, 1:17, 18:34)
  expect_true(is_efast(fit_sim))
})
