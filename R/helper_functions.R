# Helper functions
#' prepareData
#' 
#' Create data.table with input data for testing purposes
#'
#' @param cell_lines data.table with cell line info
#' @param drugs data.table with drug info
#' @param conc vector of doses
#'
#' @return data.table with input data for testing
#' 
#' @examples
#' 
#' prepareData(create_synthetic_cell_lines(), create_synthetic_drugs())
#' 
#' @export
prepareData <- function(cell_lines, drugs, conc = 10 ^ (seq(-3, 1, 0.5))) {
  df_layout <- data.table::as.data.table(merge.data.frame(cell_lines, drugs, by = NULL))
  df_layout <- add_data_replicates(df_layout)
  add_concentration(df_layout, conc)
}

#' prepareMergedData
#' 
#' Create data.table with input data containing noise for testing purposes 
#'
#' @param cell_lines data.table with cell line info
#' @param drugs data.table with drug info
#' @param noise number indicating level of noise
#'
#' @return data.table with input data for testing
#' 
#' @examples
#' 
#' prepareMergedData(create_synthetic_cell_lines(), create_synthetic_drugs())
#' 
#' @export
prepareMergedData <- function(cell_lines, drugs, noise = 0.1) {
  df <- prepareData(cell_lines, drugs)
  generate_response_data(df, noise)
}

#' prepareComboMergedData
#' 
#' Create data.table with input combination data containing noise for testing purposes 
#'
#' @param cell_lines data.table with cell line info
#' @param drugs data.table with drug info
#' @param drugsIdx1 numeric vector of ids for primary drug
#' @param drugsIdx2 numeric vector of ids for secondary drug
#' @param concentration numeric vector of doses
#' @param noise number indicating level of noise
#' @param modifyDf2 Boolean indicating if the table should me modified to keep reverse
#' single agent data
#'
#' @return data.table with input data for testing
#' 
#' @examples
#' 
#' prepareComboMergedData(create_synthetic_cell_lines(), create_synthetic_drugs())
#' 
#' @export
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

#' prepareCodilutionData
#' 
#' Create data.table with input codilution data containing noise for testing purposes 
#'
#' @param df data.table object with experiment data
#' @param df_layout data.table object with experiment design containing drugs and cell line data
#' @return data.table with input data for testing
#' 
#' @examples
#' 
#' cell_lines <- create_synthetic_cell_lines()
#' drugs <- create_synthetic_drugs()
#' df_layout <- prepareData(cell_lines[seq_len(2), ], drugs[seq_len(4), ])
#' df <- cbind(drugs[1, , drop = FALSE],
#'             df_layout[, "Concentration", drop = FALSE])
#' prepareCodilutionData(df, df_layout)
#' 
#' @export
prepareCodilutionData <- function(df, df_layout) {
  colnames(df) <- paste0(colnames(df), "_2")
  df_2 <- cbind(df_layout, df)
  df_2 <- df_2[df_2$DrugName != df_2$DrugName_2, ]
  rows <- which(df_2$Concentration_2 > 0)
  cols <- c("Concentration", "Concentration_2")
  df_2[rows, (cols) := lapply(.SD, function(x) x / 2), .SDcols = cols]
  df_2
}
