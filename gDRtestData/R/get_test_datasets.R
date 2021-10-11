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
#' @return named vector of absolute paths
#' @export
#'
#' @author Kamil Foltyński \email{kamil.foltynski@@contractors.roche.com}
#'
#' @rdname get_test_dataset_paths
get_test_dataset_paths <- function() {
  datasets_dir <- system.file("testdata", package = "gDRtestData", mustWork = TRUE)
  config <- yaml::read_yaml(get_config_path())
  out <- vapply(config, function(x) {
    checkmate::assert_string(x$file, min.chars = 1)
    fpath <- file.path(datasets_dir, x$file)
    checkmate::assert_file_exists(fpath)
    fpath
  }, FUN.VALUE = character(1))
  names(out) <- sapply(config, function(x) x$name)
  out
}


