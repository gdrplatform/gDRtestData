#' @export
generate_small_no_noise <- function() {
  cell_lines <- create_synthetic_cell_lines()
  drugs <- create_synthetic_drugs()
  e_inf <- generate_e_inf(drugs, cell_lines)
  ec50 <- generate_ec50(drugs, cell_lines)
  hill_coef <- generate_hill_coef(drugs, cell_lines)

  df_layout <- merge(cell_lines[2:11,], drugs[2:11,], by = NULL)
  df_layout <- add_data_replicates(df_layout)
  df_layout <- add_concentration(df_layout)

  df_merged_data <- generate_response_data(df_layout, 0)

  mae <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)
  mae
}
