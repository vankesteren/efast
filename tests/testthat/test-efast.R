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
  fit_efa <- efa_esem(data = lavaan::HolzingerSwineford1939[,7:15], M = 3)
  expect_true(is_efa(fit_efa))
  expect_true(inherits(efast_loadings(fit_efa), "loadings"))
})


test_that("EFAST fitting on big data works", {
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

})

test_that("Data generation works and efast thereon works", {
  set.seed(45)
  simdat <- simulate_efast()
  fit_sim <- efast_hemi(simdat, M = 4, 1:17, 18:34)
  expect_true(is_efast(fit_sim))
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
