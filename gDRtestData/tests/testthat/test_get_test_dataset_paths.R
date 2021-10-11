testthat::test_that("config exists", {
  config_path <- get_config_path()
  testthat::expect_true(tools::file_ext(config_path) %in% c("yml", "yaml"))
  testthat::expect_true(file.exists(config_path))
})


testthat::test_that("", {
  testthat::expect_type(get_test_dataset_paths(), "character")
  testthat::expect_gt(length(get_test_dataset_paths()), 0L)
})