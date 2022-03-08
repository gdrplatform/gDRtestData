testthat::test_that("test dataset paths exist", {
  testthat::expect_type(get_test_dataset_paths(), "character")
  testthat::expect_gt(length(get_test_dataset_paths()), 0L)
  dir.create(tmp_dir <- tempfile())
  on.exit(unlink(tmp_dir, TRUE, TRUE), add = TRUE)
  testthat::expect_error(get_test_dataset_paths(datasets_dir = ""),
                         "Assertion on 'datasets_dir' failed: Must have at least 1 characters.")
  testthat::expect_error(get_test_dataset_paths(datasets_dir = character(0)),
                         "Assertion on 'datasets_dir' failed: Must have length 1.")
  testthat::expect_error(get_test_dataset_paths(datasets_dir = tmp_dir, pattern = 2),
                         "Assertion on 'pattern' failed: Must be of type 'string', not 'double'.")
})
