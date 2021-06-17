#' Add data replicates
#'
#' @param df_layout data.frame that should contains the cell line,
#' drug, concentration, and replicate columns along with the annotations that needs to be propagated
#'
#' @return data.frame with replicates
#' @export
#'
add_data_replicates <- function(df_layout) {
  df_layout_duplicates <- Reduce(rbind, list(df_layout)[rep(1, times = 3)])
  barcode <- do.call(c, lapply(paste0("plate_", 1:3), function(x) rep(x, nrow(df_layout))))
  cbind(Barcode = barcode, df_layout_duplicates)
}

#' Add concentrations
#'
#' @param df_layout data.frame that should contains the cell line,
#' drug, concentration, and replicate columns along with the annotations that needs to be propagated
#' @param concentrations vector of numeric concentrations that will be added to df_layout
#'
#' @return data.frame with concentrations
#' @export
#'
add_concentration <- function(df_layout, concentrations = 10 ^ (seq(-3, 1, .5))) {
  df_layout <- merge(df_layout, data.frame(Concentration = c(0, 0, concentrations)), by = NULL)
  df_layout
}



#' Generate response data
#'
#' @param df_layout data.frame that should contains the cell line,
#' drug, concentration, and replicate columns along with the annotations that needs to be propagated
#' @param noise_level numeric scalar with the level of noise added to the data
#'
#' @return data.frame with response data
#' @export
#'
generate_response_data <- function(df_layout, noise_level = .1) {
  hill_coef <- generate_hill_coef(create_synthetic_drugs(), create_synthetic_cell_lines())
  ec50 <- generate_ec50(create_synthetic_drugs(), create_synthetic_cell_lines())
  e_inf <- generate_e_inf(create_synthetic_drugs(), create_synthetic_cell_lines())

  set.seed(2)
  df_layout$ReadoutValue <- round(100 * pmax(
    apply(df_layout, 1, function(x)
      e_inf[x["Gnumber"], x["clid"]] + (1 - e_inf[x["Gnumber"], x["clid"]]) *
        (ec50[x["Gnumber"], x["clid"]] ^ hill_coef[x["Gnumber"], x["clid"]] /
           (as.numeric(x["Concentration"]) ^ hill_coef[x["Gnumber"], x["clid"]] +
              ec50[x["Gnumber"],x["clid"]] ^ hill_coef[x["Gnumber"], x["clid"]]))) +
      (noise_level * runif(nrow(df_layout)) - (noise_level / 2)),  # add some noise
    0.01 * runif(nrow(df_layout)) + .005), # avoid hard 0 values
    1)
  df_layout$BackgroundValue <- 0
  df_layout$Duration <- 72

  #   df_layout$Gnumber <- factor(df_layout$Gnumber)
  #   df_layout$DrugName <- factor(df_layout$DrugName)
  #   df_layout$drug_moa <- factor(df_layout$drug_moa)
  #   levels(df_layout$Gnumber) <- c(levels(df_layout$Gnumber), "vehicle")
  df_layout$Gnumber[df_layout$Concentration == 0] <- "vehicle"
  #   levels(df_layout$DrugName) <- c(levels(df_layout$DrugName), "vehicle")
  df_layout$DrugName[df_layout$Concentration == 0] <- "vehicle"
  #   levels(df_layout$drug_moa) <- c(levels(df_layout$drug_moa), "vehicle")
  df_layout$drug_moa[df_layout$Concentration == 0] <- "vehicle"


  if ("Gnumber_2" %in% colnames(df_layout)) { # combo data
    df_layout$ReadoutValue <- df_layout$ReadoutValue *
      apply( df_layout, 1, function(x)
        e_inf[x["Gnumber_2"], x["clid"]] + (1 - e_inf[x["Gnumber_2"], x["clid"]])*
          (ec50[x["Gnumber_2"], x["clid"]] ^ hill_coef[x["Gnumber_2"], x["clid"]] /
             (as.numeric(x["Concentration_2"]) ^ hill_coef[x["Gnumber_2"], x["clid"]] +
                ec50[x["Gnumber_2"],x["clid"]] ^ hill_coef[x["Gnumber_2"], x["clid"]])))


    # df_layout$Gnumber_2 <- factor(df_layout$Gnumber_2)
    # df_layout$DrugName_2 <- factor(df_layout$DrugName_2)
    # df_layout$drug_moa_2 <- factor(df_layout$drug_moa_2)

    # levels(df_layout$Gnumber_2) <- c(levels(df_layout$Gnumber_2), "vehicle")
    df_layout$Gnumber_2[df_layout$Concentration_2 == 0] <- "vehicle"
    # levels(df_layout$DrugName_2) <- c(levels(df_layout$DrugName_2), "vehicle")
    df_layout$DrugName_2[df_layout$Concentration_2 == 0] <- "vehicle"
    # levels(df_layout$drug_moa_2) <- c(levels(df_layout$drug_moa_2), "vehicle")
    df_layout$drug_moa_2[df_layout$Concentration_2 == 0] <- "vehicle"
  }


  if ("Gnumber_3" %in% colnames(df_layout)) { # combo data
    df_layout$ReadoutValue <- df_layout$ReadoutValue *
      apply( df_layout, 1, function(x)
        e_inf[x["Gnumber_3"], x["clid"]] + (1 - e_inf[x["Gnumber_3"], x["clid"]]) *
          (ec50[x["Gnumber_3"], x["clid"]] ^ hill_coef[x["Gnumber_3"], x["clid"]] /
             (as.numeric(x["Concentration_3"]) ^ hill_coef[x["Gnumber_3"], x["clid"]] +
                ec50[x["Gnumber_3"], x["clid"]] ^ hill_coef[x["Gnumber_3"], x["clid"]])))

    df_layout$Gnumber_3[df_layout$Concentration_3 == 0] <- "vehicle"
    df_layout$DrugName_3[df_layout$Concentration_3 == 0] <- "vehicle"
    df_layout$drug_moa_3[df_layout$Concentration_3 == 0] <- "vehicle"
  }

  df_layout
}


#' Add data with day 0
#'
#' @param df_merged_data data.frame with merged data
#' @param noise_level numeric scalar with the level of noise added to the data
#'
#' @return
#' @export
#'
add_day0_data <- function(df_merged_data, noise_level = .05) {
  set.seed(2)
  df_Day0 = unique(df_merged_data[df_merged_data$Concentration == 0 &
                                    ifelse(array("Concentration_2", nrow(df_merged_data)) %in% colnames(df_merged_data),
                                           df_merged_data$Concentration_2 == 0, TRUE),])

  df_Day0$ReadoutValue <- df_Day0$ReadoutValue / 2 ^ (df_Day0$Duration / df_Day0$ReferenceDivisionTime)
  df_Day0$ReadoutValue <- round(df_Day0$ReadoutValue * (1 - noise_level / 2 + noise_level * runif(nrow(df_Day0))), 1)

  df_Day0$Duration <- 0
  df_Day0$Barcode <- "plate_0"

  df_merged_data <- rbind(df_merged_data, df_Day0)
  df_merged_data
}
