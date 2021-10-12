#' get_config_path
#' 
#' Return path to config in \code{yml} format containing list of test datasets
#'
#' @return absolute path to config file
#' @export
#'
#' @rdname get_test_datasets
get_config_path <- function() {
  system.file("testdata", "synthetic_list.yml", package = "gDRtestData", mustWork = TRUE)
}

#' get_test_dataset_paths
#' 
#' Returns named vector of absolute paths to test datasets.
#'
#' @param datasets_dir path to directory with datasets (default \code{NULL}).
#' If \code{NULL}, then \code{inst/testdata} directory from \code{gDRtestData} will be used.
#' @param config_path path to config i \code{yml} format (default \code{NULL}).
#' If \code{NULL}, then config from \link{get_config_path} will be used.
#'
#' @return named vector of absolute paths
#' @export
#'
#' @author Kamil Foltyński \email{kamil.foltynski@@contractors.roche.com}
#'
#' @rdname get_test_dataset_paths
get_test_dataset_paths <- function(datasets_dir = NULL, config_path = NULL) {
  if (is.null(datasets_dir))
    datasets_dir <- system.file("testdata", package = "gDRtestData", mustWork = TRUE)
  checkmate::assert_string(datasets_dir, min.chars = 1)
  checkmate::assert_directory_exists(datasets_dir)
  if (is.null(config_path))
    config_path <- get_config_path()
  checkmate::assert_string(config_path, min.chars = 1)
  checkmate::assert_file_exists(config_path)
  config <- yaml::read_yaml(config_path)
  out <- vapply(config, function(x) {
    checkmate::assert_string(x$file, min.chars = 1)
    fpath <- file.path(datasets_dir, x$file)
    checkmate::assert_file_exists(fpath)
    fpath
  }, FUN.VALUE = character(1))
  names(out) <- vapply(config, function(x) x$name, FUN.VALUE = character(1), USE.NAMES = FALSE)
  out
}


