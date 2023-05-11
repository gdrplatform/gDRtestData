#' Add data replicates
#'
#' @param df_layout data.table that should contains the cell line,
#' drug, concentration, and replicate columns along with the annotations that needs to be propagated
#'
#'
#'
#' @examples
#' 
#' cell_lines <- create_synthetic_cell_lines()
#' add_data_replicates(cell_lines)
#' 
#' @return data.table with replicates
#' @export
#'
add_data_replicates <- function(df_layout) {
  df_layout_duplicates <- Reduce(rbind, list(df_layout)[rep(1, times = 3)])
  barcode <- do.call(c, lapply(paste0("plate_", seq_len(3)), function(x) rep(x, nrow(df_layout))))
  cbind(Barcode = barcode, df_layout_duplicates)
}


#' Add concentrations
#'
#' @param df_layout data.table that should contains the cell line,
#' drug, concentration, and replicate columns along with the annotations that needs to be propagated
#' @param concentrations vector of numeric concentrations that will be added to df_layout
#'
#' @examples
#' 
#' cell_lines <- create_synthetic_cell_lines()
#' add_concentration(cell_lines)
#' 
#' @return data.table with concentrations
#' @export
#'
add_concentration <- function(df_layout, concentrations = 10 ^ (seq(-3, 1, 0.5))) {
  df_layout <- data.table::setDT(
    merge.data.frame(df_layout,
                     data.table::data.table(Concentration = c(0, 0, concentrations)),
                     by = NULL))
  df_layout
}


#' Generate response data
#'
#' @param df_layout data.table that should contains the cell line,
#' drug, concentration, and replicate columns along with the annotations that needs to be propagated
#' @param noise_level numeric scalar with the level of noise added to the data
#'
#'
#' @examples
#' 
#' cell_lines <- create_synthetic_cell_lines()
#' drugs <- create_synthetic_drugs()
#' gDRtestData:::prepareData(cell_lines[seq_len(2), ], drugs[seq_len(4), ])
#' 
#'
#' @return data.table with response data
#' @export
#'
generate_response_data <- function(df_layout, noise_level = 0.1) {
  
  drugs <- create_synthetic_drugs()
  cell_lines <- create_synthetic_cell_lines()
  hill_coef <- generate_hill_coef(drugs, cell_lines)
  ec50 <- generate_ec50(drugs, cell_lines)
  e_inf <- generate_e_inf(drugs, cell_lines)
  
  df_layout$ReadoutValue <- round(100 * pmax(
    getReadoutCoef(df_layout, e_inf, ec50, hill_coef) +
      (noise_level * stats::runif(nrow(df_layout)) - (noise_level / 2)),  # add some noise
    0.01 * stats::runif(nrow(df_layout)) + 0.005), # avoid hard 0 values
    1)
  df_layout$BackgroundValue <- 0
  df_layout$Duration <- 72
  df_layout <- introduceVehicle(df_layout)
  
  if ("Gnumber_2" %in% colnames(df_layout)) { # combo data
    df_layout <- introduceGNum(df_layout, e_inf, ec50, hill_coef, "_2")
  }
  
  if ("Gnumber_3" %in% colnames(df_layout)) { # combo data
    df_layout <- introduceGNum(df_layout, e_inf, ec50, hill_coef, "_3")
  }
  
  df_layout
}

getReadoutCoef <- function(df, e_inf, ec50, hill_coef, suffix = "") {
  apply(df, 1, function(x) {
    clid <- x["clid"]
    gnum <- x[paste0("Gnumber", suffix)]
    
    e_inf_val <- e_inf[gnum, clid]
    ec50_val <- ec50[gnum, clid]
    hill_val <- hill_coef[gnum, clid]
    concentration <- as.numeric(x["Concentration"])
    
    e_inf_val + (1 - e_inf_val) * (ec50_val ^ hill_val / (concentration ^ hill_val + ec50_val ^ hill_val)) 
  })
}

introduceVehicle <- function(df, suffix = "") {
  zeroIdx <- df[, paste0("Concentration", suffix)] == 0
  
  df[zeroIdx, paste0("Gnumber", suffix)] <- "vehicle"
  df[zeroIdx, paste0("DrugName", suffix)] <- "vehicle"
  df[zeroIdx, paste0("drug_moa", suffix)] <- "vehicle"
  
  df
}

introduceGNum <- function(df, e_inf, ec50, hill_coef, suffix) {
  df$ReadoutValue <- df$ReadoutValue * getReadoutCoef(df, e_inf, ec50, hill_coef, suffix)
  df <- introduceVehicle(df, suffix)
  
  df
}

#' Add data with day 0
#'
#' @param df_merged data.table with merged data
#' @param noise_level numeric scalar with the level of noise added to the data
#'
#' @examples
#' 
#' cell_lines <- create_synthetic_cell_lines()
#' drugs <- create_synthetic_drugs()
#' data <- gDRtestData:::prepareData(cell_lines[seq_len(2), ], drugs[seq_len(4), ])
#' data$Duration <- 72
#' data$ReadoutValue <- 0
#' add_day0_data(data)
#'
#' @return data.table with day0 data
#' @export
#'
add_day0_data <- function(df_merged, noise_level = 0.05) {
  cond <- ifelse(
    array("Concentration_2", nrow(df_merged)) %in% colnames(df_merged),
    df_merged$Concentration_2 == 0, 
    TRUE
  )
  df_Day0 <- unique(df_merged[df_merged$Concentration == 0 & cond, ])
  
  df_Day0$ReadoutValue <- df_Day0$ReadoutValue / 2 ^ (df_Day0$Duration / df_Day0$ReferenceDivisionTime)
  coef <- (1 - noise_level / 2 + noise_level * stats::runif(nrow(df_Day0)))
  df_Day0$ReadoutValue <- round(df_Day0$ReadoutValue * coef, 1)
  
  df_Day0$Duration <- 0
  df_Day0$Barcode <- "plate_0"
  
  df_merged <- rbind(df_merged, df_Day0)
  df_merged
}
