#### AUXILIARY FUNCTIONS ###

se_from_testdata <-
  function(key_valuesL,
           testDataNumb = 1,
           subdirName = "data",
           log_file = tempfile()) {

    if (!requireNamespace("gDR", quietly = TRUE)) {
      stop("Package 'gDR' needed for this function to work. Please install it.",
           call. = FALSE)
    }

    testDataDir <-
      system.file(package = "gDRtestData", "testdata", paste0(subdirName, testDataNumb))
    lRef <- gDR::read_ref_data(testDataDir)
    normSE = gDR::normalize_SE(lRef$df_raw_data, log_file, key_values = key_valuesL[[testDataNumb]])
    avgSE = gDR::average_SE(normSE)
    metricsSE = gDR::metrics_SE(avgSE)
    metricsSE
  }

#### MAIN ####
rootDir <- system.file(package = "gDRtestData", "testdata")
seL_file <- file.path(rootDir, "seL.rds")
log_file <- file.path(rootDir, "log.txt")

key_valuesL <-
  list(NULL,
       list(E2 = 1e-4),
       list(E2 = 1e-4),
       NULL,
       NULL,
       NULL,
       NULL,
       NULL,
       NULL,
       NULL,
       NULL,
       NULL,
       list(E2 = 1e-4))


# new datasets
seL <- lapply(1:13, function(x) {
  print(x)
  se_from_testdata(key_valuesL, testDataNumb = x)
})
seL <- lapply(seq_len(length(seL)), function(x){
  attr(seL[[x]], "synthetic") <- FALSE
  if (x %in% c(8,9,11,12)){
    attr(seL[[x]], "synthetic") <- TRUE
  }
  seL[[x]]
})
saveRDS(seL, seL_file)
