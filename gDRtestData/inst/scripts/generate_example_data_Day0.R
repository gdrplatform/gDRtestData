
# only for testing purpuses not displayed as example in the visualization



library(SummarizedExperiment)
devtools::load_all('../../../../BumpyMatrix')
devtools::load_all('../../../../gDRutils/gDRutils')
devtools::load_all('../../../../gDRwrapper/gDRwrapper')
devtools::load_all('../../../../gDRcore/gDRcore')
library(reshape2)

source('functions_generate_data.R')



#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the 1st test set: no noise
df_layout <- merge(CellLines[2:11,], Drugs[2:11,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout, 0)

finalSE_1_no_noise <- process_data_to_SE(df_merged_data)


df_merged_data_day0 <- add_day0_data(df_merged_data, 0)
finalSE_1_no_noise_day0 <- process_data_to_SE(df_merged_data_day0)


# testing
df_divTime = data.frame(calculated_div_time = colMeans(assay(finalSE_1_no_noise_day0,6)))
df_divTime = merge(colData(finalSE_1_no_noise_day0), df_divTime, by = 0)

all(abs(df_divTime$ReferenceDivisionTime - df_divTime$calculated_div_time)<.06)

dt = convert_se_assay_to_dt(finalSE_1_no_noise,'Averaged')
dt_day0 = convert_se_assay_to_dt(finalSE_1_no_noise_day0,'Averaged')
dt_values = merge(dt_day0, dt, by = c('rId', 'cId', 'Concentration'))

all(dt_values$RelativeViability.y - dt_values$RelativeViability.x == 0)
all(abs(dt_values$GRvalue.y - dt_values$GRvalue.x)<2e-3)


dt_test <- test_accuracy(finalSE_1_no_noise_day0, e_inf, ec50, hill_coef)
print(dt_test)
# test: 
print(apply(abs(dt_test) < c(1e-3, 2.2e-3, 0.04, 0.015, 1e-4), 1, all))

write.table(df_merged_data_day0, '../testdata/synthdata_small_wDay0_no_noise_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_1_no_noise, '../testdata/finalSE_small_wDay0_no_noise.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the 1st test set with ligand as reference
df_layout <- merge(CellLines[2:6,], Drugs[2:5,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout, 0)
df_merged_data$Ligand <- 0.1
df_merged_data2 <- df_merged_data[df_merged_data$Gnumber %in% c('vehicle', 'G00002', 'G00003'),]
df_merged_data2$Ligand <- 0
df_merged_data2$ReadoutValue <- 105 - pmax(0, pmin(104, (105-df_merged_data2$ReadoutValue) ** 1.1))
df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0('CL000', 11:12)] <- 
    0.8 * df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0('CL000', 11:12)]
df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0('CL000', 13:14)] <- 
    0.5 * df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0('CL000', 13:14)]
df_merged_data2$ReadoutValue <- round(df_merged_data2$ReadoutValue,1)
df_merged_data2$Barcode <- paste0(df_merged_data2$Barcode, '1')
df_merged_data <- rbind(df_merged_data, df_merged_data2)

df_merged_data_day0 <- add_day0_data(df_merged_data, 0)
df_merged_data_day0 = df_merged_data_day0[df_merged_data_day0$Ligand == 0.1 | df_merged_data_day0$Duration>0,]

finalSE_1_Ligand <- process_data_to_SE(df_merged_data, override_untrt_controls = c(Ligand = 0.1))
finalSE_1_Ligand_day0 <- process_data_to_SE(df_merged_data_day0, override_untrt_controls = c(Ligand = 0.1))

dt_test <- test_accuracy(finalSE_1_Ligand_day0[rowData(finalSE_1_Ligand_day0)$Ligand > 0,], e_inf, ec50, hill_coef)
print(dt_test)
# test: 
print(apply(abs(dt_test) < c(1e-3, 3e-3, 0.031, 0.015, 1e-4), 1, all))

dt = convert_se_assay_to_dt(finalSE_1_Ligand,'Averaged')
dt_day0 = convert_se_assay_to_dt(finalSE_1_Ligand_day0,'Averaged')
dt_values = merge(dt_day0, dt, by = c('rId', 'cId', 'Concentration'))

all(dt_values$RelativeViability.y - dt_values$RelativeViability.x == 0)
all(abs(dt_values$GRvalue.y - dt_values$GRvalue.x)<2e-3)

write.table(df_merged_data, '../testdata/synthdata_wLigand_wDay0_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_1_Ligand, '../testdata/finalSE_wLigand_wDay0.RDS', compress = FALSE)

