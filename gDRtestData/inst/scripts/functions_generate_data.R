
#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generic parameters for 60 drugs across 15 cell lines (3 tissues)
CellLines <- data.frame(
    clid = paste0('CL000', 11:25),
    CellLineName = paste0('cellline_', LETTERS[1:15]),
    Tissue = sort(paste0('tissue_', array(letters[c(24, 24:26, 26)], 15))),
    ReferenceDivisionTime = seq(22,80,4)
)
CellLines <- rbind(CellLines, CellLines, CellLines, CellLines, CellLines, CellLines)
CellLines$clid <- paste0('CL000', 9+(seq_len(nrow(CellLines))))
CellLines$CellLineName <- paste0(CellLines$CellLineName, sort(array(LETTERS[1:15], 90)))
CellLines$Tissue[16:40] <- 'tissue_w'
CellLines$Tissue[41:50] <- 'tissue_v'

Drugs <- data.frame(
    Gnumber = paste0('G00', 11:50),
    DrugName = paste0('drug_', 11:50),
    drug_moa = sort(paste0('moa_', array(LETTERS[c(1, 1:6, 6)], 40)))
)
Drugs <- rbind(Drugs, Drugs, Drugs, Drugs, Drugs, Drugs)
Drugs$Gnumber <- sprintf('G00%03i', seq_len(nrow(Drugs)))
Drugs$DrugName <- sprintf('drug_%03i', seq_len(nrow(Drugs)))
Drugs$drug_moa[-1:-80] <- sort(paste0('moa_', array(LETTERS[1:24], 160)))


set.seed(2)
hill_coef <- matrix(1.8 + runif(nrow(Drugs) * nrow(CellLines)), nrow(Drugs), nrow(CellLines))
colnames(hill_coef) <- CellLines$clid
rownames(hill_coef) <- Drugs$Gnumber

set.seed(2)
ec50 <- matrix( runif(nrow(Drugs) * nrow(CellLines))-.5, nrow(Drugs), nrow(CellLines)) + 
          matrix(sort(rep(seq(-1.2,0,.3),8)),nrow(Drugs),nrow(CellLines)) +
          t(matrix(seq(-.4,0,.1), nrow(CellLines), nrow(Drugs)))

ec50[, CellLines$Tissue=='tissue_x'] <- ec50[, CellLines$Tissue=='tissue_x'] +
          -.5 -.3*runif(sum(CellLines$Tissue=='tissue_x'))
ec50[Drugs$drug_moa %in% c('moa_A', 'moa_B'), CellLines$Tissue=='tissue_y'] <- 
    ec50[Drugs$drug_moa %in% c('moa_A', 'moa_B'), CellLines$Tissue=='tissue_y'] +
          -.4*runif(sum(CellLines$Tissue=='tissue_y')) -.4
ec50[Drugs$drug_moa %in% c('moa_C', 'moa_D'), CellLines$Tissue=='tissue_z'] <- 
    ec50[Drugs$drug_moa %in% c('moa_C', 'moa_D'), CellLines$Tissue=='tissue_z'] +
              .5*runif(sum(CellLines$Tissue=='tissue_z')) - 1
ec50[Drugs$drug_moa %in% 'moa_E',] <- ec50[Drugs$drug_moa %in% 'moa_E',] + .6 + .5*runif(sum(Drugs$drug_moa %in% 'moa_E'))
ec50 <- 10 **(ec50)

colnames(ec50) <- CellLines$clid
rownames(ec50) <- Drugs$Gnumber
quantile(ec50)

set.seed(2)
e_inf <- matrix( .5*runif(nrow(Drugs) * nrow(CellLines)), nrow(Drugs), nrow(CellLines))  +
          t(matrix(seq(0,.2,.05), nrow(CellLines), nrow(Drugs)))  

e_inf[, CellLines$Tissue=='tissue_x'] <- e_inf[, CellLines$Tissue=='tissue_x'] +
          .3*runif(sum(CellLines$Tissue=='tissue_x')) + .1
e_inf[Drugs$drug_moa %in% c('moa_A', 'moa_C'), CellLines$Tissue=='tissue_y'] <- 
    e_inf[Drugs$drug_moa %in% c('moa_A', 'moa_C'), CellLines$Tissue=='tissue_y'] +
          -.2*runif(sum(CellLines$Tissue=='tissue_y')) -.2
e_inf[Drugs$drug_moa %in% c('moa_C', 'moa_E'), CellLines$Tissue=='tissue_z'] <- 
    e_inf[Drugs$drug_moa %in% c('moa_C', 'moa_E'), CellLines$Tissue=='tissue_z'] -
              .5*runif(sum(CellLines$Tissue=='tissue_z'))
e_inf[Drugs$drug_moa %in% c('moa_B', 'moa_E'), CellLines$Tissue=='tissue_x'] <- 
    e_inf[Drugs$drug_moa %in% c('moa_B', 'moa_E'), CellLines$Tissue=='tissue_x'] -
              .5*runif(sum(CellLines$Tissue=='tissue_x'))
e_inf[Drugs$drug_moa %in% 'moa_F',] <- e_inf[Drugs$drug_moa %in% 'moa_F',] + .3 + .2*runif(sum(Drugs$drug_moa %in% 'moa_F'))

e_inf <- matrix( pmin(.89, pmax(.01, e_inf)) + runif(nrow(Drugs) * nrow(CellLines))*.1, nrow(Drugs), nrow(CellLines))
colnames(e_inf) <- CellLines$clid
rownames(e_inf) <- Drugs$Gnumber
quantile(e_inf)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 




#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# function to generate the data
# df_layout should contains the cell line, drug, concentration, and replicate columns 
#     along with the annotations that needs to be propagated

add_data_replicates <- function(df_layout) {
  df_layout <- rbind(cbind(Barcode = 'plate_1', df_layout),
              cbind(Barcode = 'plate_2', df_layout),
              cbind(Barcode = 'plate_3', df_layout))
}

add_concentration <- function(df_layout, Concentrations = 10**(seq(-3,1,.5))) {
  df_layout <- merge(df_layout, data.frame(Concentration = c(0, 0, Concentrations)), by = NULL)
}



generate_response_data <- function(df_layout, noise_level = .1) {
  set.seed(2)
  df_layout$ReadoutValue <- round(100 * pmax(
        apply( df_layout, 1, function(x)
          e_inf[x['Gnumber'],x['clid']] + (1-e_inf[x['Gnumber'],x['clid']])*
              (ec50[x['Gnumber'],x['clid']] ** hill_coef[x['Gnumber'],x['clid']] / 
                (as.numeric(x['Concentration']) ** hill_coef[x['Gnumber'],x['clid']] +
                  ec50[x['Gnumber'],x['clid']] ** hill_coef[x['Gnumber'],x['clid']]))) +
          (noise_level*runif(nrow(df_layout)) - (noise_level/2)),  # add some noise
        0.01*runif(nrow(df_layout))+.005), # avoid hard 0 values
      1)
  quantile(df_layout$ReadoutValue)
  df_layout$BackgroundValue <- 0
  df_layout$Duration <- 72

#   df_layout$Gnumber <- factor(df_layout$Gnumber)
#   df_layout$DrugName <- factor(df_layout$DrugName)
#   df_layout$drug_moa <- factor(df_layout$drug_moa)
#   levels(df_layout$Gnumber) <- c(levels(df_layout$Gnumber), 'vehicle')
  df_layout$Gnumber[df_layout$Concentration %in% 0] <- 'vehicle'
#   levels(df_layout$DrugName) <- c(levels(df_layout$DrugName), 'vehicle')
  df_layout$DrugName[df_layout$Concentration %in% 0] <- 'vehicle'
#   levels(df_layout$drug_moa) <- c(levels(df_layout$drug_moa), 'vehicle')
  df_layout$drug_moa[df_layout$Concentration %in% 0] <- 'vehicle'
  

  if ('Gnumber_2' %in% colnames(df_layout)) { # combo data
    df_layout$ReadoutValue <- df_layout$ReadoutValue *
      apply( df_layout, 1, function(x)
        e_inf[x['Gnumber_2'],x['clid']] + (1-e_inf[x['Gnumber_2'],x['clid']])*
            (ec50[x['Gnumber_2'],x['clid']] ** hill_coef[x['Gnumber_2'],x['clid']] / 
              (as.numeric(x['Concentration_2']) ** hill_coef[x['Gnumber_2'],x['clid']] +
                ec50[x['Gnumber_2'],x['clid']] ** hill_coef[x['Gnumber_2'],x['clid']])))


    # df_layout$Gnumber_2 <- factor(df_layout$Gnumber_2)
    # df_layout$DrugName_2 <- factor(df_layout$DrugName_2)
    # df_layout$drug_moa_2 <- factor(df_layout$drug_moa_2)
    
    # levels(df_layout$Gnumber_2) <- c(levels(df_layout$Gnumber_2), 'vehicle')
    df_layout$Gnumber_2[df_layout$Concentration_2 == 0] <- 'vehicle'
    # levels(df_layout$DrugName_2) <- c(levels(df_layout$DrugName_2), 'vehicle')
    df_layout$DrugName_2[df_layout$Concentration_2 == 0] <- 'vehicle'
    # levels(df_layout$drug_moa_2) <- c(levels(df_layout$drug_moa_2), 'vehicle')
    df_layout$drug_moa_2[df_layout$Concentration_2 %in% 0] <- 'vehicle'
  }


  if ('Gnumber_3' %in% colnames(df_layout)) { # combo data
    df_layout$ReadoutValue <- df_layout$ReadoutValue *
      apply( df_layout, 1, function(x)
        e_inf[x['Gnumber_3'],x['clid']] + (1-e_inf[x['Gnumber_3'],x['clid']])*
            (ec50[x['Gnumber_3'],x['clid']] ** hill_coef[x['Gnumber_3'],x['clid']] / 
              (as.numeric(x['Concentration_3']) ** hill_coef[x['Gnumber_3'],x['clid']] +
                ec50[x['Gnumber_3'],x['clid']] ** hill_coef[x['Gnumber_3'],x['clid']])))

    df_layout$Gnumber_3[df_layout$Concentration_3 == 0] <- 'vehicle'
    df_layout$DrugName_3[df_layout$Concentration_3 == 0] <- 'vehicle'
    df_layout$drug_moa_3[df_layout$Concentration_3 %in% 0] <- 'vehicle'
  }

  return(df_layout)
}


add_day0_data <- function(df_merged_data, noise_level = .05) {
  set.seed(2)
  df_Day0 = unique(df_merged_data[df_merged_data$Concentration == 0 &
    ifelse(array('Concentration_2',nrow(df_merged_data)) %in% colnames(df_merged_data), 
        df_merged_data$Concentration_2 == 0, T),])

  df_Day0$ReadoutValue = df_Day0$ReadoutValue / 2**(df_Day0$Duration / df_Day0$ReferenceDivisionTime)
  df_Day0$ReadoutValue = round(df_Day0$ReadoutValue * (1 - noise_level/2 + noise_level*runif(nrow(df_Day0))),1)

  df_Day0$Duration = 0
  df_Day0$Barcode = 'plate_0'

  df_merged_data = rbind(df_merged_data, df_Day0)
}

process_data_to_SE2 <- function(df_data, override_untrt_controls = NULL) {
  se <- gDR::create_SE2(df_data, override_untrt_controls = override_untrt_controls)
  normSE <- gDR::normalize_SE2(se)  
  avgSE <- gDR::average_SE2(normSE)
  metricsSE <- gDR::fit_SE2(avgSE)
  finalSE <- gDR::add_codrug_group_SE(metricsSE)
}

# function to test accuracy of the fitted metrics based on the model
test_accuracy <- function(finalSE, e_inf, ec50, hill_coef) {
  
  # by passing until 'flatten' in debugged

  # dt <- gDRutils::convert_se_assay_to_dt(finalSE, 'Metrics')
  # gDRutils::flatten(dt, groups = c('normalization_type', 'fit_source'), 
  #       wide_cols = gDRutils::get_header('response_metrics'))

  # df_QC <- rbind(quantile(acast(dt, Gnumber ~ clid, value.var = 'E_inf') - 
  #     e_inf[ rowData(finalSE)$Gnumber, colData(finalSE)$clid ], c(.05, .5, .95)),
  #   quantile(log10(acast(dt, Gnumber ~ clid, value.var = 'EC50')) - 
  #     log10(ec50[ rowData(finalSE)$Gnumber, colData(finalSE)$clid ]), c(.05, .5, .95)),
  #   quantile(acast(dt, Gnumber ~ clid, value.var = 'h_RV') - 
  #     hill_coef[ rowData(finalSE)$Gnumber, colData(finalSE)$clid ], c(.05, .5, .95)),
  #   quantile( (acast(dt, Gnumber ~ clid, value.var = 'h_RV') - 
  #     hill_coef[ rowData(finalSE)$Gnumber, colData(finalSE)$clid ])[
  #       acast(dt, Gnumber ~ clid, value.var = 'EC50') < 3 & 
  #       acast(dt, Gnumber ~ clid, value.var = 'E_inf') < .8 
  #     ], c(.05, .5, .95)),
  #   1-quantile(acast(dt, Gnumber ~ clid, value.var = 'RV_r2') , c(.05, .5, .95))
  # )

  df_QC = as.data.frame( matrix(1, 5, 3)) # to remove 

  rownames(df_QC) <- c('delta_einf', 'delta_ec50', 'delta_hill', 'd_hill_fitted', '1_r2')
  return(df_QC)
}
#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


