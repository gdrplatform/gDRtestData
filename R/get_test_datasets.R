#' get_test_dataset_paths
#'
#' Returns named vector of absolute paths to test datasets.
#'
#' @param datasets_dir path to directory with datasets (default \code{NULL}).
#' If \code{NULL}, then \code{inst/testdata} directory from \code{gDRtestData} will be used.
#' @param pattern used to: (1) filter to qs2 files from the dataset_dir path
#'        and (2) prettify the labels of the files
#' @keywords generate_test_data
#'
#' @return named vector of absolute paths
#'
#' @examples
#'
#' get_test_dataset_paths()
#' path <- system.file("testdata", package = "gDRtestData", mustWork = TRUE)
#' get_test_dataset_paths(path)
#'
#' @export
#'
#' @author Kamil Foltyński \email{kamil.foltynski@@contractors.roche.com}
#'
#' @rdname get_test_dataset_paths
get_test_dataset_paths <-
  function(datasets_dir = NULL,
           pattern = "finalMAE_") {
   if (is.null(datasets_dir)) {
     datasets_dir <-
       system.file("testdata", package = "gDRtestData", mustWork = TRUE)
   }
   checkmate::assert_string(datasets_dir, min.chars = 1)
   checkmate::assert_directory_exists(datasets_dir)

   checkmate::assert_string(pattern, min.chars = 1)

   epaths <- list.files(datasets_dir, pattern = paste0(pattern, ".*\\.qs2$"), full.names = TRUE)
   enames <- gsub(pattern, "", gsub("\\.qs2$", "", basename(epaths)))
   structure(epaths, names = enames)
 }
