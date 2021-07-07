#' Test accuracy of the fitted metrics based on the model
#'
#' @param se SummarizedExperiment for testing
#' @param e_inf e_inf
#' @param ec50 ec50
#' @param hill_coef hill_coef
#'
#' @return data.frame with quality metrics
#' @export
#'
test_accuracy <- function(se, e_inf, ec50, hill_coef) {

  dt <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  dt <- gDRutils::flatten(dt, groups = c("normalization_type", "fit_source"),
                          wide_cols = gDRutils::get_header("response_metrics"))
  colnames(dt) <- gDRutils::prettify_flat_metrics(colnames(dt), FALSE)
  df_QC <- rbind(quantile(acast(dt, Gnumber ~ clid, value.var = "E_inf") -
                            e_inf[ rowData(se)$Gnumber, colData(se)$clid ], c(.05, .5, .95)),
                 quantile(log10(acast(dt, Gnumber ~ clid, value.var = "EC50")) -
                            log10(ec50[rowData(se)$Gnumber, colData(se)$clid]), c(.05, .5, .95)),
                 quantile(acast(dt, Gnumber ~ clid, value.var = "h_RV") -
                            hill_coef[ rowData(se)$Gnumber, colData(se)$clid ], c(.05, .5, .95)),
                 quantile( (acast(dt, Gnumber ~ clid, value.var = "h_RV") -
                              hill_coef[ rowData(se)$Gnumber, colData(se)$clid ])[
                                acast(dt, Gnumber ~ clid, value.var = "EC50") < 3 &
                                  acast(dt, Gnumber ~ clid, value.var = "E_inf") < .8
                                ], c(.05, .5, .95)),
                 1 - quantile(acast(dt, Gnumber ~ clid, value.var = "RV_r2") , c(.05, .5, .95))
  )

  rownames(df_QC) <- c("delta_einf", "delta_ec50", "delta_hill", "d_hill_fitted", "1_r2")
  df_QC
}
