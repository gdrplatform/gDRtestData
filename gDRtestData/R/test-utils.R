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
  rows <- unique(SummarizedExperiment::rowData(se)$Gnumber)
  cols <- SummarizedExperiment::colData(se)$clid
  quart <- c(.05, .5, .95)
    
  df_QC <- rbind(
    stats::quantile(acastVar(dt, "E_inf") - e_inf[rows, cols], quart),
    stats::quantile(log10(acastVar(dt, "EC50")) - log10(ec50[rows, cols]), quart),
    stats::quantile(acastVar(dt, "h_RV") - hill_coef[rows, cols], quart),
    stats::quantile((acastVar(dt, "h_RV") - hill_coef[rows, cols])[
      acastVar(dt, "EC50") < 3 & acastVar(dt, "E_inf") < .8], quart
    ),
    1 - stats::quantile(acastVar(dt, "RV_r2"), quart)
  )

  rownames(df_QC) <- c("delta_einf", "delta_ec50", "delta_hill", "d_hill_fitted", "1_r2")
  df_QC
}

acastVar <- function(dt, var) {
  reshape2::acast(dt, Gnumber ~ clid, value.var = var)
}
