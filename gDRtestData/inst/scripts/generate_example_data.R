# generating three types of example/test data
#   1. 10 drugs (3 differents gcsi_moa) by 10 lines (3 tissues); single agent - no noise in the data
#   2. 10 drugs (3 differents gcsi_moa) by 10 lines (3 tissues); single agent
#   3. 40 drugs (6 differents gcsi_moa) by 15 lines (3 tissues); single agent
#   4. 10 drugs x 2 co-treatment (1 drugs at 2 doses) by 10 lines (3 tissues);
#   4. 4 drugs x 1 drug (co-dilution) by 8 lines (3 tissues);
#   4. 6 drugs x 3 drugs (matrix co-treatments) by 4 lines (3 tissues);

library(SummarizedExperiment)
devtools::load_all('../../../../BumpyMatrix')
devtools::load_all('../../../../gDRutils/gDRutils')
devtools::load_all('../../../../gDRwrapper/gDRwrapper')
devtools::load_all('../../../../gDRcore/gDR')
library(reshape2)



#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generic parameters for 60 drugs across 15 cell lines (3 tissues)
CellLines = data.frame(
    clid = paste0('CL000', 11:25),
    CellLineName = paste0('cellline_', LETTERS[1:15]),
    Tissue = sort(paste0('tissue_', array(letters[c(24, 24:26, 26)], 15))),
    ReferenceDivisionTime = seq(22,80,4)
)
CellLines = rbind(CellLines, CellLines, CellLines, CellLines, CellLines, CellLines)
CellLines$clid = paste0('CL000', 9+(seq_len(nrow(CellLines))))
CellLines$CellLineName = paste0(CellLines$CellLineName, sort(array(LETTERS[1:15], 90)))
CellLines$Tissue[16:40] = 'tissue_w'
CellLines$Tissue[41:50] = 'tissue_v'

Drugs = data.frame(
    Gnumber = paste0('G00', 11:50),
    DrugName = paste0('drug_', 11:50),
    gcsi_moa = sort(paste0('moa_', array(LETTERS[c(1, 1:6, 6)], 40)))
)
Drugs = rbind(Drugs, Drugs, Drugs, Drugs, Drugs, Drugs)
Drugs$Gnumber = sprintf('G00%03i', seq_len(nrow(Drugs)))
Drugs$DrugName = sprintf('drug_%03i', seq_len(nrow(Drugs)))
Drugs$gcsi_moa[-1:-80] = sort(paste0('moa_', array(LETTERS[1:24], 160)))


set.seed(2)
hill_coef = matrix(1.8 + runif(nrow(Drugs) * nrow(CellLines)), nrow(Drugs), nrow(CellLines))
colnames(hill_coef) = CellLines$clid
rownames(hill_coef) = Drugs$Gnumber

set.seed(2)
ec50 <- matrix( runif(nrow(Drugs) * nrow(CellLines))-.5, nrow(Drugs), nrow(CellLines)) + 
          matrix(sort(rep(seq(-1.2,0,.3),8)),nrow(Drugs),nrow(CellLines)) +
          t(matrix(seq(-.4,0,.1), nrow(CellLines), nrow(Drugs)))

ec50[, CellLines$Tissue=='tissue_x'] <- ec50[, CellLines$Tissue=='tissue_x'] +
          -.5 -.3*runif(sum(CellLines$Tissue=='tissue_x'))
ec50[Drugs$gcsi_moa %in% c('moa_A', 'moa_B'), CellLines$Tissue=='tissue_y'] <- 
    ec50[Drugs$gcsi_moa %in% c('moa_A', 'moa_B'), CellLines$Tissue=='tissue_y'] +
          -.4*runif(sum(CellLines$Tissue=='tissue_y')) -.4
ec50[Drugs$gcsi_moa %in% c('moa_C', 'moa_D'), CellLines$Tissue=='tissue_z'] <- 
    ec50[Drugs$gcsi_moa %in% c('moa_C', 'moa_D'), CellLines$Tissue=='tissue_z'] +
              .5*runif(sum(CellLines$Tissue=='tissue_z')) - 1
ec50[Drugs$gcsi_moa %in% 'moa_E',] <- ec50[Drugs$gcsi_moa %in% 'moa_E',] + .6 + .5*runif(sum(Drugs$gcsi_moa %in% 'moa_E'))
ec50 <- 10 **(ec50)

colnames(ec50) = CellLines$clid
rownames(ec50) = Drugs$Gnumber
quantile(ec50)

set.seed(2)
e_max = matrix( .5*runif(nrow(Drugs) * nrow(CellLines)), nrow(Drugs), nrow(CellLines))  +
          t(matrix(seq(0,.2,.05), nrow(CellLines), nrow(Drugs)))  

e_max[, CellLines$Tissue=='tissue_x'] <- e_max[, CellLines$Tissue=='tissue_x'] +
          .3*runif(sum(CellLines$Tissue=='tissue_x')) + .1
e_max[Drugs$gcsi_moa %in% c('moa_A', 'moa_C'), CellLines$Tissue=='tissue_y'] <- 
    e_max[Drugs$gcsi_moa %in% c('moa_A', 'moa_C'), CellLines$Tissue=='tissue_y'] +
          -.2*runif(sum(CellLines$Tissue=='tissue_y')) -.2
e_max[Drugs$gcsi_moa %in% c('moa_C', 'moa_E'), CellLines$Tissue=='tissue_z'] <- 
    e_max[Drugs$gcsi_moa %in% c('moa_C', 'moa_E'), CellLines$Tissue=='tissue_z'] -
              .5*runif(sum(CellLines$Tissue=='tissue_z'))
e_max[Drugs$gcsi_moa %in% c('moa_B', 'moa_E'), CellLines$Tissue=='tissue_x'] <- 
    e_max[Drugs$gcsi_moa %in% c('moa_B', 'moa_E'), CellLines$Tissue=='tissue_x'] -
              .5*runif(sum(CellLines$Tissue=='tissue_x'))
e_max[Drugs$gcsi_moa %in% 'moa_F',] <- e_max[Drugs$gcsi_moa %in% 'moa_F',] + .3 + .2*runif(sum(Drugs$gcsi_moa %in% 'moa_F'))

e_max <- matrix( pmin(.89, pmax(.01, e_max)) + runif(nrow(Drugs) * nrow(CellLines))*.1, nrow(Drugs), nrow(CellLines))
colnames(e_max) = CellLines$clid
rownames(e_max) = Drugs$Gnumber
quantile(e_max)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 




#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# function to generate the data
# df_raw should contains the cell line, drug, concentration, and replicate columns 
#     along with the annotations that needs to be propagated

add_data_replicates <- function(df_raw) {
  df_raw <- rbind(cbind(Barcode = 'plate_1', df_raw),
              cbind(Barcode = 'plate_2', df_raw),
              cbind(Barcode = 'plate_3', df_raw))
}
add_concentration <- function(df_raw) {
  Concentrations = 10**(seq(-3,1,.5))
  df_raw <- merge(df_raw, data.frame(Concentration = c(0, 0, Concentrations)), by = NULL)
}

generate_response_data <- function(df_raw, noise_level = .1) {
  set.seed(2)
  df_raw$ReadoutValue <- round(100 * pmax(
        apply( df_raw, 1, function(x)
          e_max[x['Gnumber'],x['clid']] + (1-e_max[x['Gnumber'],x['clid']])*
              (ec50[x['Gnumber'],x['clid']] ** hill_coef[x['Gnumber'],x['clid']] / 
                (as.numeric(x['Concentration']) ** hill_coef[x['Gnumber'],x['clid']] +
                  ec50[x['Gnumber'],x['clid']] ** hill_coef[x['Gnumber'],x['clid']]))) +
          (noise_level*runif(nrow(df_raw)) - (noise_level/2)),  # add some noise
        0.01*runif(nrow(df_raw))+.005), # avoid hard 0 values
      1)
  quantile(df_raw$ReadoutValue)
  df_raw$BackgroundValue <- 0
  df_raw$Duration <- 72

#   df_raw$Gnumber = factor(df_raw$Gnumber)
#   df_raw$DrugName = factor(df_raw$DrugName)
#   df_raw$gcsi_moa = factor(df_raw$gcsi_moa)
#   levels(df_raw$Gnumber) <- c(levels(df_raw$Gnumber), 'vehicle')
  df_raw$Gnumber[df_raw$Concentration %in% 0] <- 'vehicle'
#   levels(df_raw$DrugName) <- c(levels(df_raw$DrugName), 'vehicle')
  df_raw$DrugName[df_raw$Concentration %in% 0] <- 'vehicle'
#   levels(df_raw$gcsi_moa) <- c(levels(df_raw$gcsi_moa), 'vehicle')
  df_raw$gcsi_moa[df_raw$Concentration %in% 0] <- 'vehicle'
  

  if ('Gnumber_2' %in% colnames(df_raw)) { # combo data
    df_raw$ReadoutValue = df_raw$ReadoutValue *
      apply( df_raw, 1, function(x)
        e_max[x['Gnumber_2'],x['clid']] + (1-e_max[x['Gnumber_2'],x['clid']])*
            (ec50[x['Gnumber_2'],x['clid']] ** hill_coef[x['Gnumber_2'],x['clid']] / 
              (as.numeric(x['Concentration_2']) ** hill_coef[x['Gnumber_2'],x['clid']] +
                ec50[x['Gnumber_2'],x['clid']] ** hill_coef[x['Gnumber_2'],x['clid']])))


    # df_raw$Gnumber_2 = factor(df_raw$Gnumber_2)
    # df_raw$DrugName_2 = factor(df_raw$DrugName_2)
    # df_raw$gcsi_moa_2 = factor(df_raw$gcsi_moa_2)
    
    # levels(df_raw$Gnumber_2) <- c(levels(df_raw$Gnumber_2), 'vehicle')
    df_raw$Gnumber_2[df_raw$Concentration_2 == 0] <- 'vehicle'
    # levels(df_raw$DrugName_2) <- c(levels(df_raw$DrugName_2), 'vehicle')
    df_raw$DrugName_2[df_raw$Concentration_2 == 0] <- 'vehicle'
    # levels(df_raw$gcsi_moa_2) <- c(levels(df_raw$gcsi_moa_2), 'vehicle')
    df_raw$gcsi_moa_2[df_raw$Concentration %in% 0] <- 'vehicle'
  }

  return(df_raw)
}


process_data_to_SE2 <- function(df_raw, key_values = NULL) {
  se <- gDR::create_SE2(df_raw, key_values = key_values)
  normSE <- gDR::normalize_SE2(se)  
  avgSE <- gDR::average_SE2(normSE)
  metricsSE <- gDR::fit_SE2(avgSE)
#   finalSE <- gDR::add_codrug_group_SE(metricsSE)
}

process_data_to_SE <- function(df_raw, key_values = NULL) {
  normSE <- gDR::normalize_SE(df_raw, key_values)  
  avgSE <- gDR::average_SE(normSE)
  metricsSE <- gDR::metrics_SE(avgSE)
  finalSE <- gDR::add_codrug_group_SE(metricsSE)
}
# function to test accuracy of the fitted metrics based on the model
test_accuracy <- function(finalSE) {
  dt = assay_to_dt(finalSE, 'Metrics', merge_metrics = TRUE)

  df_QC = rbind(quantile(acast(dt, Gnumber ~ clid, value.var = 'e_max') - 
      e_max[ rowData(finalSE)$Gnumber, colData(finalSE)$clid ], c(.05, .5, .95)),
    quantile(log10(acast(dt, Gnumber ~ clid, value.var = 'ec50')) - 
      log10(ec50[ rowData(finalSE)$Gnumber, colData(finalSE)$clid ]), c(.05, .5, .95)),
    quantile(acast(dt, Gnumber ~ clid, value.var = 'h_RV') - 
      hill_coef[ rowData(finalSE)$Gnumber, colData(finalSE)$clid ], c(.05, .5, .95)),
    quantile( (acast(dt, Gnumber ~ clid, value.var = 'h_RV') - 
      hill_coef[ rowData(finalSE)$Gnumber, colData(finalSE)$clid ])[
        acast(dt, Gnumber ~ clid, value.var = 'ec50') < 3 & 
        acast(dt, Gnumber ~ clid, value.var = 'e_max') < .8 
      ], c(.05, .5, .95)),
    1-quantile(acast(dt, Gnumber ~ clid, value.var = 'RV_r2') , c(.05, .5, .95))
  )
  rownames(df_QC) = c('delta_emax', 'delta_ec50', 'delta_hill', 'd_hill_fitted', '1_r2')
  return(df_QC)
}
#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 





#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the 1st test set: no noise
#   only for testing purpuses not displayed as example
df_raw = merge(CellLines[2:11,], Drugs[2:11,], by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_merged_data = generate_response_data(df_raw, 0)

# finalSE_1_no_noise = process_data_to_SE(df_merged_data)
finalSE_1_no_noise = process_data_to_SE2(df_merged_data)
dt_test = test_accuracy(finalSE_1_no_noise)
print(dt_test)
# test: 
apply(abs(dt_test) < c(1e-3, 2e-3, 0.02, 0.015, 1e-4),1,all)

saveRDS(finalSE_1_no_noise, '../testdata/finalSE_small_no_noise.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the 1st test set with noise
df_raw = merge(CellLines[2:11,], Drugs[2:11,], by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_merged_data = generate_response_data(df_raw)
finalSE_1 = process_data_to_SE(df_merged_data)

dt_test = test_accuracy(finalSE_1)
print(dt_test)
# test: 
apply(abs(dt_test) < c(0.5, 0.1, 1.5, 1.2, 0.05),1,all)

write.table(df_merged_data, '../testdata/synthdata_small_rawdata.tsv', quote = F, row.names = F, sep = '\t')
write.table(assay_to_dt(finalSE_1, 3)[,-1:-2], '../testdata/synthdata_small_avg.tsv', quote = F, row.names = F, sep = '\t')
write.table(assay_to_dt(finalSE_1, 5, T)[,-1:-2], '../testdata/synthdata_small_metrics.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_1, '../testdata/finalSE_small.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the 1st test set with noise and ligand as reference
df_raw = merge(CellLines[2:8,], Drugs[2:8,], by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_merged_data = generate_response_data(df_raw)
df_merged_data$E2 = 0.1
df_merged_data2 = df_merged_data[df_merged_data$Gnumber %in% c('vehicle', 'G00002', 'G00003', 'G00004'),]
df_merged_data2$E2 = 0
df_merged_data2$ReadoutValue = 105 - pmax(0, pmin(104, (105-df_merged_data2$ReadoutValue) ** 1.1))
df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0('CL000', 11:14)] = 
    0.8 * df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0('CL000', 11:14)]
df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0('CL000', 15:16)] = 
    0.5 * df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0('CL000', 15:16)]
df_merged_data2$ReadoutValue = round(df_merged_data2$ReadoutValue,1)
df_merged_data2$Barcode = paste0(df_merged_data2$Barcode, '1')
df_merged_data = rbind(df_merged_data, df_merged_data2)

finalSE_1 = process_data_to_SE(df_merged_data)
# finalSE_1 = process_data_to_SE(df_merged_data, key_values = c(E2 = 0.1))

write.table(df_merged_data, '../testdata/synthdata_E2_ref_rawdata.tsv', quote = F, row.names = F, sep = '\t')
write.table(assay_to_dt(finalSE_1, 3)[,-1:-2], '../testdata/synthdata_E2_ref_avg.tsv', quote = F, row.names = F, sep = '\t')
write.table(assay_to_dt(finalSE_1, 5, T)[,-1:-2], '../testdata/synthdata_E2_ref_metrics.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_1, '../testdata/finalSE_1.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the 2nd (medium size) test set with single agent
df_raw = merge(CellLines[1:15,], Drugs[1:40,], by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_merged_data = generate_response_data(df_raw)
finalSE_2 = process_data_to_SE(df_merged_data)

dt_test = test_accuracy(finalSE_2)
print(dt_test)
# test: 
apply(abs(dt_test) < c(0.5, 0.2, 2.5, 1.2, 0.3),1,all)

write.table(df_merged_data, '../testdata/synthdata_medium_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_2, '../testdata/finalSE_medium.RDS', compress = FALSE)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the 3rd (many lines) test set with single agent
df_raw = merge(CellLines, Drugs[1:40,], by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_merged_data = generate_response_data(df_raw)
finalSE_2 = process_data_to_SE(df_merged_data)

dt_test = test_accuracy(finalSE_2)
print(dt_test)
# test: 
apply(abs(dt_test) < c(0.5, 0.2, 2.5, 1.2, 0.3),1,all)

write.table(df_merged_data, '../testdata/synthdata_many_lines_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_2, '../testdata/finalSE_many_lines.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the test set with single agent (many drugs)
df_raw = merge(CellLines[1:10,], Drugs, by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_merged_data = generate_response_data(df_raw)
finalSE_2 = process_data_to_SE(df_merged_data)

dt_test = test_accuracy(finalSE_2)
print(dt_test)
# test: 
apply(abs(dt_test) < c(0.5, 0.2, 2.5, 1.2, 0.3),1,all)

write.table(df_merged_data, '../testdata/synthdata_many_drugs_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_2, '../testdata/finalSE_many_drugs.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the 3rd test set with combo (two single dose)
df_raw = merge(CellLines[2:4,], Drugs[2:4,], by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_2 = cbind(Drugs[c(1,1,1),], Concentration = c(0, .2, 1))
colnames(df_2) = paste0(colnames(df_2), '_2')
df_raw_2 = merge(df_raw, df_2, by = NULL)

df_merged_data = generate_response_data(df_raw_2)

gDR::create_SE2(df_merged_data)

finalSE_combo = process_data_to_SE2(df_merged_data)
# finalSE_combo = process_data_to_SE(df_merged_data)

saveRDS(finalSE_combo, '../testdata/finalSE_combo_2dose.RDS', compress = FALSE)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the test set with combo (unique dose; many drug)
df_raw = merge(CellLines[2:4,], Drugs, by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_2 = cbind(Drugs[c(1,1),], Concentration = c(0, 2))
colnames(df_2) = paste0(colnames(df_2), '_2')
df_raw_2 = merge(df_raw, df_2, by = NULL)

df_merged_data = generate_response_data(df_raw_2)
finalSE_combo = process_data_to_SE(df_merged_data)

saveRDS(finalSE_combo, '../testdata/finalSE_combo_1dose_many_drugs.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the 4th test set with combo (co-dilution)
df_raw = merge(CellLines[seq(1,15,2),], Drugs[1:4,], by = NULL)
df_raw = add_data_replicates(df_raw)
df_raw = add_concentration(df_raw)

df_2 = cbind(Drugs[c(1,1),], df_raw[,'Concentration',drop=F])
colnames(df_2) = paste0(colnames(df_2), '_2')
df_raw_2 = cbind(df_raw, df_2)
df_raw_2[df_raw_2$Concentration_2 > 0 , c('Concentration', 'Concentration_2')] = 
  df_raw_2[df_raw_2$Concentration_2 > 0 , c('Concentration', 'Concentration_2')] / 2

df_merged_data = generate_response_data(df_raw_2)
finalSE_codilution = process_data_to_SE(df_merged_data)
# add calculation for co-dilution

saveRDS(finalSE_codilution, '../testdata/finalSE_combo_codilution.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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

df_merged_data = generate_response_data(df_raw_2)
finalSE_matrix = process_data_to_SE(df_merged_data)
# add calculation for combo matrix

saveRDS(finalSE_matrix, '../testdata/finalSE_combo_matrix.RDS', compress = FALSE)