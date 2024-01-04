# Helper functions
#' @keywords internal
prepareData <- function(cell_lines, drugs, conc = 10 ^ (seq(-3, 1, 0.5))) {
  df_layout <- drugs[, as.list(cell_lines), names(drugs)]
  df_layout <- add_data_replicates(df_layout)
  add_concentration(df_layout, conc)
}

#' @keywords internal
prepareMergedData <- function(cell_lines, drugs, noise = 0.1) {
  df <- prepareData(cell_lines, drugs)
  generate_response_data(df, noise)
}

#' @keywords internal
prepareComboMergedData <- function(cell_lines, 
                                   drugs, 
                                   drugsIdx1 = 2:4,
                                   drugsIdx2 = c(26, 26, 26), 
                                   concentration = c(0, .2, 1), 
                                   noise = 0.1, 
                                   modifyDf2 = FALSE) {
  df_layout <- prepareData(cell_lines, drugs[drugsIdx1, ])
  
  df_2 <- cbind(drugs[drugsIdx2, ], Concentration = concentration)
  colnames(df_2) <- paste0(colnames(df_2), "_2")
  
  df_layout_2 <- data.table::as.data.table(merge.data.frame(df_layout, df_2, by = NULL))
  if (modifyDf2) {
    df_layout_2 <- df_layout_2[!(df_layout_2$Concentration == 0 & df_layout_2$Concentration_2 > 0), ]
  }
  
  generate_response_data(df_layout_2, noise)
}

#' @keywords internal
prepareCodilutionData <- function(df, df_layout) {
  colnames(df) <- paste0(colnames(df), "_2")
  df_2 <- cbind(df_layout, df)
  df_2 <- df_2[df_2$DrugName != df_2$DrugName_2, ]
  rows <- which(df_2$Concentration_2 > 0)
  cols <- c("Concentration", "Concentration_2")
  df_2[rows, (cols) := lapply(.SD, function(x) x / 2), .SDcols = cols]
  df_2
}
