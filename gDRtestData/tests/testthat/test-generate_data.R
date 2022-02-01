library(testthat); library(gDRtestData)

test_that("generate_small_no_noise works as expected", {
  expect_true(is(generate_small_no_noise(), "MultiAssayExperiment"))
})
