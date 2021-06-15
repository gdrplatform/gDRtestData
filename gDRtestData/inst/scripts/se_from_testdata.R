#### AUXILIARY FUNCTIONS ###

se_from_testdata <-
  function(key_valuesL,
           discard_keys_l,
           testDataNumb = 1,
           subdirName = "data",
           log_file = tempfile()) {

    if (!requireNamespace("gDRcore", quietly = TRUE)) {
      stop("Package 'gDRcore' needed for this function to work. Please install it.",
           call. = FALSE)
    }

    testDataDir <-
      system.file(package = "gDRtestData", "testdata", paste0(subdirName, testDataNumb))
    lRef <- gDRcore::read_ref_data(testDataDir)
    se <– gDRcore::create_SE(lRef$df_raw_data,
                             nested_keys = key_valuesL[[testDataNumb]],
                             override_untrt_controls = discard_keys_l[[testDataNumb]])
    normSE <- gDRcore::normalize_SE(se)
    avgSE <- gDRcore::average_SE(normSE)
    metricsSE <- gDRcore::fit_SE(avgSE)
    finalmetricsSE <- gDRcore::add_codrug_group_SE(metricsSE)
    finalmetricsSE
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

discard_keys_l <-
  list(
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    "Replicate",
    "Replicate",
    "Replicate",
    NULL,
    "Replicate",
    "Replicate",
    NULL
  )

# new datasets
seL <- lapply(1:13, function(x) {
  print(x)
  se_from_testdata(key_valuesL, discard_keys_l, testDataNumb = x)
})
seL <- lapply(seq_len(length(seL)), function(x){
  print(x)
  attr(seL[[x]], "synthetic") <- FALSE
  if (x %in% c(8,9,11,12)){
    attr(seL[[x]], "synthetic") <- TRUE
  }
  seL[[x]]
})
saveRDS(seL, seL_file)
