#' @export
generateNoNoiseRawData <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # generate the data for the 1st test set: no noise
  #   only for testing purpuses not displayed as example
  df_merged <- prepareMergedData(cell_lines[2:11, ], drugs[2:11, ], 0)
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "synthdata_small_no_noise_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_small_no_noise.RDS"
    )
  }
  
  invisible(mae)
}

#' @export
generateNoiseRawData <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  # generate the data for the 1st test set with noise
  df_merged <- prepareMergedData(cell_lines[2:11, ], drugs[2:11, ])
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "synthdata_small_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_small.RDS"
    )
  }
  
  invisible(mae)
}

#' @export
generateLigandData <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  # generate the data for the 1st test set with ligand as reference
  df_merged <- prepareMergedData(cell_lines[2:6, ], drugs[2:5, ], 0)
  df_merged$Ligand <- 0.1
  df_merged2 <- df_merged[df_merged$Gnumber %in% c("vehicle", "G00002", "G00003"), ]
  df_merged2$Ligand <- 0
  
  idx1 <- df_merged2$clid %in% paste0("CL000", 11:12)
  idx2 <- df_merged2$clid %in% paste0("CL000", 13:14)
  df_merged2$ReadoutValue <- 105 - pmax(0, pmin(104, (105 - df_merged2$ReadoutValue) ^ 1.1))
  df_merged2$ReadoutValue[idx1] <- 0.8 * df_merged2$ReadoutValue[idx1]
  df_merged2$ReadoutValue[idx2] <- 0.5 * df_merged2$ReadoutValue[idx2]
  df_merged2$ReadoutValue <- round(df_merged2$ReadoutValue, 1)
  
  df_merged2$Barcode <- paste0(df_merged2$Barcode, "1")
  df_merged <- rbind(df_merged, df_merged2)
  
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged, 
    override_untrt_controls = c(Ligand = 0.1),
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "finalSE_wLigand_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_wLigand.RDS"
    )
  }
  
  invisible(mae)
}

#' @export
generateMediumData <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  # generate the data for the 2nd (medium size) test set with single agent
  df_merged <- prepareMergedData(cell_lines[1:15, ], drugs[1:40, ])
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )

  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "synthdata_medium_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_medium.RDS"
    )
  }

  invisible(mae)
}

#' @export
generateManyLinesData <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  # generate the data for the 2nd (medium size) test set with single agent
  df_merged <- prepareMergedData(cell_lines, drugs[1:40, ])
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )

  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "synthdata_many_lines_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_many_lines.RDS"
    )
  }
  
  invisible(mae)
}

#' @export
generateManyDrugsData <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  # generate the data for the test set with single agent (many drugs)
  df_merged <- prepareMergedData(cell_lines[1:10, ], drugs[1:40, ])
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "synthdata_many_drugs_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_many_drugs.RDS"
    )
  }
  
  invisible(mae)
}

#' @export
generateComboNoNoiseData <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  # generate the data for the test set with combo (two single dose)
  #   co-treatment drug is only as DrugName_2
  df_merged <- prepareComboMergedData(cell_lines[2:4, ], drugs, noise = 0)
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "synthdata_combo_2dose_nonoise_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_combo_2dose_nonoise.RDS"
    )
  }
  
  invisible(mae)
}

#' @export
generateComboNoNoiseData2 <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  # generate the data for the test set with combo (two single dose)
  #   co-treatment drug is also as single agent as DrugName
  df_merged <- prepareComboMergedData(cell_lines[2:4, ], drugs, drugsIdx1 = c(2:4, 26), noise = 0)
  df_merged <- df_merged[!(df_merged$Gnumber %in% c("vehicle", drugs$Gnumber[26]) &
    df_merged$Gnumber_2 == drugs$Gnumber[26]), ]
  
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "synthdata_combo_2dose_nonoise2_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_combo_2dose_nonoise2.RDS"
    )
  }
  
  invisible(mae)
}

#' @export
generateComboNoNoiseData3 <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  # generate the data for the 3rd test set with combo (two single dose)
  #   co-treatment drug does NOT have single agent response
  df_merged <- prepareComboMergedData(
    cell_lines = cell_lines[2:4, ], 
    drugs = drugs, 
    drugsIdx1 = 2:4, 
    noise = 0, 
    modifyDf2 = TRUE
  )
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )

  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "synthdata_combo_2dose_nonoise3_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_combo_2dose_nonoise3.RDS"
    )
  }
  
  invisible(mae)
}

#' @export
generateComboManyDrugs <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  # generate the data for the test set with combo (unique dose; many drug)
  df_merged <- prepareComboMergedData(
    cell_lines = cell_lines[2:4, ], 
    drugs = drugs,
    drugsIdx1 = -1,
    drugsIdx2 = c(1, 1),
    concentration = c(0, 2)
  )
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )

  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "synthdata_combo_1dose_many_drugs_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_combo_1dose_many_drugs.RDS"
    )
  }
  
  invisible(mae)
}

#' @export
generateComboMatrixSmall <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  # generate the data with combo matrix (small, no noise)
  concentration <- 10^ (seq(-3, .5, .5))
  df_layout <- prepareData(cell_lines[7:8, ], drugs[c(4:6), ], concentration)
  df_2 <- prepareData(cell_lines[cell_lines$clid %in% df_layout$clid, ], drugs[c(21, 26), ], concentration)
  df_2 <- changeColNames(df_2, drugs, "_2")
  df_layout_2 <- merge(df_layout, df_2)
  
  df_merged <- generate_response_data(df_layout_2, 0)
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "synthdata_combo_matrix_small_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_combo_matrix_small.RDS"
    )
  }
  
  invisible(mae)
}

#' @export
generateComboMatrix <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  # generate the data with combo matrix (mid-size)
  df_layout <- prepareData(cell_lines[seq(1, 30, 4), ], drugs[c(1, 2, 11, 12, 16, 17), ])
  df_2 <- prepareData(cell_lines[cell_lines$clid %in% df_layout$clid, ], drugs[c(21, 26, 31), ])
  df_2 <- changeColNames(df_2, drugs, "_2")
  df_layout_2 <- merge(df_layout, df_2)
  
  df_merged <- generate_response_data(df_layout_2)
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "synthdata_combo_matrix_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_combo_matrix.RDS"
    )
  }
  
  invisible(mae)
}

#' @export
generateTripleComboMatrix <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  # generate the data with triple combo  (no noise)
  concentration <- 10^ (seq(-3, .5, .5))
  df_layout <- prepareData(cell_lines[7:8, ], drugs[c(4:6), ], concentration)
  
  df_2 <- prepareData(
    cell_lines[cell_lines$clid %in% df_layout$clid, ], 
    drugs[c(21, 26), ],
    c(0, concentration)
  )
  df_2 <- changeColNames(df_2, drugs, "_2")
  df_layout_2 <- merge(df_layout, df_2)
  
  df_3 <- prepareData(
    cell_lines[cell_lines$clid %in% df_layout$clid, ], 
    drugs[10, ],
    c(0, .1, 1)
  )
  df_3 <- changeColNames(df_3, drugs, "_3")
  df_layout_3 <- merge(merge(df_layout, df_2), df_3)
  
  df_merged <- generate_response_data(df_layout_3, 0)
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "synthdata_combo_triple_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_combo_triple.RDS"
    )
  }
  
  invisible(mae)
}

#' @export
generateCodilutionSmall <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  # generate the data with combo co-dilution (small)
  df_layout <- prepareData(cell_lines[1:2, ], drugs[1:4, ])

  df_2 <- cbind(drugs[1, , drop = FALSE], df_layout[, "Concentration", drop = FALSE])
  df_layout_2 <- prepareCodilutionData(df_2, df_layout)
  
  df_merged <- generate_response_data(df_layout_2, 0)
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "synthdata_combo_codilution_small_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_combo_codilution_small.RDS"
    )
  }
  
  invisible(mae)
}

#' @export
generateCodilution <- function(cell_lines, drugs, e_inf, ec50, hill_coef, save = TRUE) {
  # generate the data for the test set with combo (co-dilution)
  df_layout <- prepareData(cell_lines[seq(1, 15, 2), ], drugs[1:12, ])
  
  df_2 <- cbind(drugs[c(1, 1), ], df_layout[, "Concentration", drop = FALSE])
  df_layout_2 <- prepareCodilutionData(df_2, df_layout)

  df_merged <- generate_response_data(df_layout_2)
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    saveArtifacts(
      tsvObj = df_merged,
      tsvName = "synthdata_combo_codilution_rawdata.tsv",
      rdsObj = mae,
      rdsName = "finalMAE_combo_codilution.RDS"
    )
  }
  
  invisible(mae)
}
