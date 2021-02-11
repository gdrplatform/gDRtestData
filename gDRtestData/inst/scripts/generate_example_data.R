# generating three types of example/test data
#   1. 10 drugs (3 differents MOA) by 10 lines (3 tissues); single agent - no noise in the data
#   2. 10 drugs (3 differents MOA) by 10 lines (3 tissues); single agent
#   3. 40 drugs (6 differents MOA) by 15 lines (3 tissues); single agent
#   4. 10 drugs x 2 co-treatment (1 drugs at 2 doses) by 10 lines (3 tissues);
#   4. 4 drugs x 1 drug (co-dilution) by 8 lines (3 tissues);
#   4. 6 drugs x 3 drugs (matrix co-treatments) by 4 lines (3 tissues);

library(SummarizedExperiment)
devtools::load_all('../../../../BumpyMatrix')
devtools::load_all('../../../../gDRutils/gDRutils')
devtools::load_all('../../../../gDRwrapper/gDRwrapper')
devtools::load_all('../../../../gDRcore/gDR')
library(reshape2)

# generic functions and parameters for 60 drugs across 15 cell lines (3 tissues)
CellLines = data.frame(
    clid = paste0('CL000', 11:25),
    CellLineName = paste0('cellline_', LETTERS[1:15]),
    Tissue = sort(paste0('tissue_', array(letters[c(24, 24:26, 26)], 15))),
    ReferenceDivisionTime = seq(22,80,4)
)

Drugs = data.frame(
    Gnumber = paste0('G00', 11:50),
    DrugName = paste0('drug_', 11:50),
    MOA = sort(paste0('moa_', array(LETTERS[c(1, 1:6, 6)], 40)))
)

set.seed(2)
hill_coef = matrix(1.8 + runif(600), 40, 15)
colnames(hill_coef) = CellLines$clid
rownames(hill_coef) = Drugs$Gnumber

set.seed(2)
ec50 <- matrix( runif(600)-.5, nrow(Drugs), nrow(CellLines)) + 
          matrix(sort(rep(seq(-1.2,0,.3),8)),nrow(Drugs),nrow(CellLines)) +
          t(matrix(seq(-.4,0,.1), nrow(CellLines), nrow(Drugs)))

ec50[, CellLines$Tissue=='tissue_x'] <- ec50[, CellLines$Tissue=='tissue_x'] +
          -.5 -.3*runif(sum(CellLines$Tissue=='tissue_x'))
ec50[Drugs$MOA %in% c('moa_A', 'moa_B'), CellLines$Tissue=='tissue_y'] <- 
    ec50[Drugs$MOA %in% c('moa_A', 'moa_B'), CellLines$Tissue=='tissue_y'] +
          -.4*runif(sum(CellLines$Tissue=='tissue_y')) -.4
ec50[Drugs$MOA %in% c('moa_C', 'moa_D'), CellLines$Tissue=='tissue_z'] <- 
    ec50[Drugs$MOA %in% c('moa_C', 'moa_D'), CellLines$Tissue=='tissue_z'] +
              .5*runif(sum(CellLines$Tissue=='tissue_z')) - 1
ec50[Drugs$MOA %in% 'moa_E',] <- ec50[Drugs$MOA %in% 'moa_E',] + .6 + .5*runif(sum(Drugs$MOA %in% 'moa_E'))
ec50 <- 10 **(ec50)

colnames(ec50) = CellLines$clid
rownames(ec50) = Drugs$Gnumber
quantile(ec50)

set.seed(2)
e_max = matrix( .5*runif(600), nrow(Drugs), nrow(CellLines))  +
          t(matrix(seq(0,.2,.05), nrow(CellLines), nrow(Drugs)))  

e_max[, CellLines$Tissue=='tissue_x'] <- e_max[, CellLines$Tissue=='tissue_x'] +
          .3*runif(sum(CellLines$Tissue=='tissue_x')) + .1
e_max[Drugs$MOA %in% c('moa_A', 'moa_C'), CellLines$Tissue=='tissue_y'] <- 
    e_max[Drugs$MOA %in% c('moa_A', 'moa_C'), CellLines$Tissue=='tissue_y'] +
          -.2*runif(sum(CellLines$Tissue=='tissue_y')) -.2
e_max[Drugs$MOA %in% c('moa_C', 'moa_E'), CellLines$Tissue=='tissue_z'] <- 
    e_max[Drugs$MOA %in% c('moa_C', 'moa_E'), CellLines$Tissue=='tissue_z'] -
              .5*runif(sum(CellLines$Tissue=='tissue_z'))
e_max[Drugs$MOA %in% c('moa_B', 'moa_E'), CellLines$Tissue=='tissue_x'] <- 
    e_max[Drugs$MOA %in% c('moa_B', 'moa_E'), CellLines$Tissue=='tissue_x'] -
              .5*runif(sum(CellLines$Tissue=='tissue_x'))
e_max[Drugs$MOA %in% 'moa_F',] <- e_max[Drugs$MOA %in% 'moa_F',] + .3 + .2*runif(sum(Drugs$MOA %in% 'moa_F'))

e_max <- matrix( pmin(.89, pmax(.01, e_max)) + runif(600)*.1, nrow(Drugs), nrow(CellLines))
colnames(e_max) = CellLines$clid
rownames(e_max) = Drugs$Gnumber
quantile(e_max)





add_data_replicates <- function(df_raw) {
  df_raw <- rbind(cbind(Barcode = 'plate_1', df_raw),
              cbind(Barcode = 'plate_2', df_raw),
              cbind(Barcode = 'plate_3', df_raw))
}
add_concentration <- function(df_raw) {
  Concentrations = 10**(seq(-3,1,.5))
  df_raw <- merge(df_raw, data.frame(Concentration = c(0, 0, Concentrations)), by = NULL)
}

# function to generate the data
# df_raw should contains the cell line, drug, concentration, and replicate columns 
#     along with the annotations that needs to be propagated
generate_response_data <- function(df_raw, noise_level = .1) {
  set.seed(2)
  df_raw$ReadoutValue <- round(100 * pmax(
        apply( df_raw, 1, function(x)
          e_max[x['Gnumber'],x['clid']] + (1-e_max[x['Gnumber'],x['clid']])*
              (ec50[x['Gnumber'],x['clid']]**2 / (as.numeric(x['Concentration'])**2 +
                  ec50[x['Gnumber'],x['clid']]**2))) +
          (noise_level*runif(nrow(df_raw)) - (noise_level/2)),  # add some noise
        0.01*runif(nrow(df_raw))+.005), # avoid hard 0 values
      1)
  quantile(df_raw$ReadoutValue)
  df_raw$BackgroundValue <- 0
  df_raw$Duration <- 72

  if ('Gnumber_2' %in% colnames(df_raw)) { # combo data
    df_raw$ReadoutValue = df_raw$ReadoutValue *
      apply( df_raw, 1, function(x)
        e_max[x['Gnumber_2'],x['clid']] + (1-e_max[x['Gnumber_2'],x['clid']])*
            (ec50[x['Gnumber_2'],x['clid']]**2 / (as.numeric(x['Concentration_2'])**2 +
                ec50[x['Gnumber_2'],x['clid']]**2)))

    levels(df_raw$Gnumber_2) <- c(levels(df_raw$Gnumber_2), 'vehicle')
    df_raw$Gnumber_2[df_raw$Concentration_2 == 0] <- 'vehicle'
    levels(df_raw$DrugName_2) <- c(levels(df_raw$DrugName_2), 'vehicle')
    df_raw$DrugName_2[df_raw$Concentration_2 == 0] <- 'vehicle'
  }

  levels(df_raw$Gnumber) <- c(levels(df_raw$Gnumber), 'vehicle')
  df_raw$Gnumber[df_raw$Concentration == 0] <- 'vehicle'
  levels(df_raw$DrugName) <- c(levels(df_raw$DrugName), 'vehicle')
  df_raw$DrugName[df_raw$Concentration == 0] <- 'vehicle'
  
  return(df_raw)
}


process_data_to_SE <- function(df_raw) {
  normSE <- gDR::normalize_SE(df_raw)  
  avgSE <- gDR::average_SE(normSE)
  metricsSE <- gDR::metrics_SE(avgSE)
  finalSE <- gDR::add_codrug_group_SE(metricsSE)
}

# function to test accuracy of the fitted metrics based on the model
test_accuracy <- function(finalSE) {
  dt = assay_to_dt(finalSE, 5, merge_metrics = TRUE)

  df_QC = rbind(quantile(acast(dt, Gnumber ~ clid, value.var = 'e_max') - 
      e_max[ rowData(finalSE)$Gnumber, colData(finalSE)$clid ], c(.05, .5, .95)),
    quantile(acast(dt, Gnumber ~ clid, value.var = 'ec50') - 
      ec50[ rowData(finalSE)$Gnumber, colData(finalSE)$clid ], c(.05, .5, .95)),
    1-quantile(acast(dt, Gnumber ~ clid, value.var = 'RV_r2') , c(.05, .5, .95))
  )
  rownames(df_QC) = c('delta_emax', 'delta_ec50', '1_r2')
  return(df_QC)
}


# generate the data for the 1st test set: no noise
#   only for testing purpuses not displayed as example
df_raw = merge(CellLines[2:11,], Drugs[2:11,], by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_data = generate_response_data(df_raw, 0)
finalSE_1_no_noise = process_data_to_SE(df_data)
print(test_accuracy(finalSE_1_no_noise))
# test: 
all(abs(test_accuracy(finalSE_1_no_noise))<1e-3)

saveRDS(finalSE_1_no_noise, '../testdata/finalSE_1_no_noise.RDS', compress = FALSE)


# generate the data for the 1st test set with noise
df_raw = merge(CellLines[2:11,], Drugs[2:11,], by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_data = generate_response_data(df_raw)
finalSE_1 = process_data_to_SE(df_data)
print(test_accuracy(finalSE_1))

saveRDS(finalSE_1, '../testdata/finalSE_1.RDS', compress = FALSE)


# generate the data for the 2nd (large) test set with single agent
df_raw = merge(CellLines, Drugs, by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_data = generate_response_data(df_raw)
finalSE_2 = process_data_to_SE(df_data)
print(test_accuracy(finalSE_2))

saveRDS(finalSE_2, '../testdata/finalSE_2.RDS', compress = FALSE)


# generate the data for the 3rd test set with combo (single dose)
df_raw = merge(CellLines[2:11,], Drugs[2:11,], by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_2 = cbind(Drugs[c(1,1,1),], Concentration = c(0, .2, 1))
colnames(df_2) = paste0(colnames(df_2), '_2')
df_raw_2 = merge(df_raw, df_2, by = NULL)

df_data = generate_response_data(df_raw_2)
finalSE_combo = process_data_to_SE(df_data)

saveRDS(finalSE_combo, '../testdata/finalSE_combo.RDS', compress = FALSE)


# generate the data for the 4th test set with combo (co-dilution)
df_raw = merge(CellLines[seq(1,15,2),], Drugs[1:4,], by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_2 = cbind(Drugs[c(1,1),], df_raw[,'Concentration',drop=F])
colnames(df_2) = paste0(colnames(df_2), '_2')
df_raw_2 = cbind(df_raw, df_2)
df_raw_2[df_raw_2$Concentration_2 > 0 , c('Concentration', 'Concentration_2')] = 
  df_raw_2[df_raw_2$Concentration_2 > 0 , c('Concentration', 'Concentration_2')] / 2

df_data = generate_response_data(df_raw_2)
finalSE_codilution = process_data_to_SE(df_data)
# add calculation for co-dilution

saveRDS(finalSE_codilution, '../testdata/finalSE_codilution.RDS', compress = FALSE)


# generate the data for the 5th test set with combo (matrix)
df_raw = merge(CellLines[seq(1,15,4),], Drugs[c(1,2,11,12,16,17),], by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_2 = merge(CellLines[CellLines$clid %in% df_raw$clid,], Drugs[c(21,26,31),], by = NULL)
df_2 = add_data_replicates(df_2)
df_2 = add_concentration(df_2)
colnames(df_2)[colnames(df_2) %in% c(colnames(Drugs),'Concentration')] = 
    paste0(colnames(df_2)[colnames(df_2) %in% c(colnames(Drugs),'Concentration')], '_2')

df_raw_2 = merge(df_raw, df_2)

df_data = generate_response_data(df_raw_2)
finalSE_matrix = process_data_to_SE(df_data)
# add calculation for combo matrix

saveRDS(finalSE_matrix, '../testdata/finalSE_matrix.RDS', compress = FALSE)