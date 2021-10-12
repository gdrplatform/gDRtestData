testthat::test_that("config exists", {
  config_path <- get_config_path()
  testthat::expect_true(tools::file_ext(config_path) %in% c("yml", "yaml"))
  testthat::expect_true(file.exists(config_path))
})


testthat::test_that("test dataset paths exist", {
  testthat::expect_type(get_test_dataset_paths(), "character")
  testthat::expect_gt(length(get_test_dataset_paths()), 0L)
  dir.create(tmp_dir <- tempfile())
  on.exit(unlink(tmp_dir, TRUE, TRUE), add = TRUE)
  testthat::expect_error(get_test_dataset_paths(datasets_dir = ""),
                         "Assertion on 'datasets_dir' failed: Must have at least 1 characters.")
  testthat::expect_error(get_test_dataset_paths(datasets_dir = character(0)),
                         "Assertion on 'datasets_dir' failed: Must have length 1.")
  testthat::expect_error(get_test_dataset_paths(datasets_dir = tmp_dir),
                         "Assertion on 'fpath' failed: File does not exist")
  testthat::expect_error(get_test_dataset_paths(config_path = tempfile()),
                         "Assertion on 'config_path' failed: File does not exist")
  testthat::expect_error(get_test_dataset_paths(config_path = ""),
                         "Assertion on 'config_path' failed: Must have at least 1 characters.")
  testthat::expect_error(get_test_dataset_paths(config_path = character(0)),
                         "Assertion on 'config_path' failed: Must have length 1.")
  dummy_config <- list(list(name = "dummy_dataset", file = "some/non/existing/dataset"))
  dummy_config_path <- tempfile()
  on.exit(unlink(dummy_config_path, TRUE, TRUE), add = TRUE)
  yaml::write_yaml(dummy_config, dummy_config_path)
  
  testthat::expect_error(get_test_dataset_paths(dirname(dummy_config_path), basename(dummy_config_path)),
                         "Assertion on 'config_path' failed: File does not exist:")
  testthat::expect_error(get_test_dataset_paths(config_path = dummy_config_path),
                         "Assertion on 'fpath' failed: File does not exist:")
})