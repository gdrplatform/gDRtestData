# generating three types of example/test data
#
#  The data have no day0 information.
#  dataset with * are to be imported as example for visualization
#
#  "small_no_noise"           10 drugs (3 differents drug_moa) by 10 lines (3 tissues); single agent - no noise in the data
#  "small"                *    10 drugs (3 differents drug_moa) by 10 lines (3 tissues); single agent
#  "wLigand"              *    3 drugs by 4 lines (3 tissues); "Ligand = 0.1" as reference; single agent
#  "medium"               *    40 drugs (6 differents drug_moa) by 15 lines (3 tissues); single agent
#  "many_lines"           *    150 drugs (6 differents drug_moa) by 10 lines (3 tissues); single agent
#  "many_drugs"           *    150 drugs (6 differents drug_moa) by 10 lines (3 tissues); single agent
#  "combo_2dose_nonoise"  *    3 drugs x 2 co-treatment (1 drug at 2 doses) by 3 lines; co-treatment drug is only as DrugName_2
#  "combo_2dose_nonoise2"      3 drugs x 2 co-treatment (1 drug at 2 doses) by 3 lines; co-treatment drug is also as single agent as DrugName
#  "combo_2dose_nonoise3"      3 drugs x 2 co-treatment (1 drug at 2 doses) by 3 lines; co-treatment drug does NOT have single agent response
#  "combo_1dose_many_drugs" *  149 drugs x 1 drug (1 dose) by 3 lines;
#  "combo_matrix_small"        3 x 2 drugs (matrix) for 2 cell lines; no noise
#  "combo_matrix"           *  6 x 3 drugs (matrix) for 8 cell lines
#  "combo_triple"           *  2 x 3 x 2 drugs (few doses) for 4 cell lines
#  "combo_codilution_small"    4 x 1 drugs (co-dilution) for 2 cell lines; no noise
#  "combo_codilution"          12 x 1 drugs (co-dilution) for 8 cell lines; no noise

library(SummarizedExperiment)
devtools::load_all('../../../../BumpyMatrix')
devtools::load_all('../../../../gDRutils/gDRutils')
devtools::load_all('../../../../gDRwrapper/gDRwrapper')
devtools::load_all('../../../../gDRcore/gDR')
library(reshape2)

source('functions_generate_data.R')



#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the 1st test set: no noise
#   only for testing purpuses not displayed as example
df_layout <- merge(CellLines[2:11,], Drugs[2:11,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout, 0)

finalSE_1_no_noise <- process_data_to_SE(df_merged_data)

# test accurarcy of the processing and fitting (no noise => low tolerance)
dt_test <- test_accuracy(finalSE_1_no_noise, e_inf, ec50, hill_coef)
print(dt_test)
# test: 
print(apply(abs(dt_test) < c(1e-3, 2.2e-3, 0.04, 0.015, 1e-4), 1, all))

write.table(df_merged_data, '../testdata/synthdata_small_no_noise_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_1_no_noise, '../testdata/finalSE_small_no_noise.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the 1st test set with noise
df_layout <- merge(CellLines[2:11,], Drugs[2:11,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout)
finalSE_1 <- process_data_to_SE(df_merged_data)

# test accurarcy of the processing and fitting (noise => medium tolerance)
dt_test <- test_accuracy(finalSE_1, e_inf, ec50, hill_coef)
print(dt_test)
# test: 
print(apply(abs(dt_test) < c(0.5, 0.1, 1.5, 1.2, 0.05), 1, all))

write.table(df_merged_data, '../testdata/synthdata_small_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_1, '../testdata/finalSE_small.RDS', compress = FALSE)


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

finalSE_1_Ligand <- process_data_to_SE(df_merged_data, override_untrt_controls = c(Ligand = 0.1))

# test accurarcy of the processing and fitting for Ligand = 0.1 (noise => medium tolerance)
dt_test <- test_accuracy(finalSE_1_Ligand[rowData(finalSE_1_Ligand)$Ligand > 0,], e_inf, ec50, hill_coef)
print(dt_test)
# test: 
print(apply(abs(dt_test) < c(1e-3, 3e-3, 0.031, 0.015, 1e-4), 1, all))
# test fit quality for Ligand = 0 and that delta(e_inf) < 0
dt_test <- test_accuracy(finalSE_1_Ligand[rowData(finalSE_1_Ligand)$Ligand == 0,], e_inf, ec50, hill_coef)
print(dt_test)
# test:
print(apply(dt_test[c('delta_einf', '1_r2'),] < c(-0.15, 1e-4), 1, all))

write.table(df_merged_data, '../testdata/finalSE_wLigand_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_1_Ligand, '../testdata/finalSE_wLigand.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the 2nd (medium size) test set with single agent
df_layout <- merge(CellLines[1:15,], Drugs[1:40,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout)
finalSE_2 <- process_data_to_SE(df_merged_data)

# test accurarcy of the processing and fitting 
dt_test <- test_accuracy(finalSE_2, e_inf, ec50, hill_coef)
print(dt_test)
# test: 
print(apply(abs(dt_test) < c(0.5, 0.2, 2.5, 1.2, 0.3), 1, all))

write.table(df_merged_data, '../testdata/synthdata_medium_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_2, '../testdata/finalSE_medium.RDS', compress = FALSE)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the 3rd (many lines) test set with single agent
df_layout <- merge(CellLines, Drugs[1:40,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout)
finalSE_3 <- process_data_to_SE(df_merged_data)

# test accurarcy of the processing and fitting 
dt_test <- test_accuracy(finalSE_3, e_inf, ec50, hill_coef)
print(dt_test)
# test: 
print(apply(abs(dt_test) < c(0.5, 0.2, 2.5, 1.2, 0.3), 1, all))

write.table(df_merged_data, '../testdata/synthdata_many_lines_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_3, '../testdata/finalSE_many_lines.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the test set with single agent (many drugs)
df_layout <- merge(CellLines[1:10,], Drugs, by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout)
finalSE_4 <- process_data_to_SE(df_merged_data)

# test accurarcy of the processing and fitting 
dt_test <- test_accuracy(finalSE_4, e_inf, ec50, hill_coef)
print(dt_test)
# test: 
print(apply(abs(dt_test) < c(0.5, 0.2, 2.5, 1.2, 0.3), 1, all))

write.table(df_merged_data, '../testdata/synthdata_many_drugs_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_4, '../testdata/finalSE_many_drugs.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the test set with combo (two single dose)
#   co-treatment drug is only as DrugName_2
df_layout <- merge(CellLines[2:4,], Drugs[2:4,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(Drugs[c(26,26,26),], Concentration = c(0, .2, 1))
colnames(df_2) <- paste0(colnames(df_2), '_2')
df_layout_2 <- merge(df_layout, df_2, by = NULL)

df_merged_data <- generate_response_data(df_layout_2, 0)
finalSE_combo <- process_data_to_SE(df_merged_data)

# test the single agent response
dt_test <- test_accuracy(finalSE_combo[rowData(finalSE_combo)$Concentration_2 == 0, ], e_inf, ec50, hill_coef)
print(dt_test)
# test: 
print(apply(abs(dt_test) < c(1e-3, 2e-3, 0.02, 0.015, 1e-4), 1, all))
# test the assignment of drug_combinations
print(length(metadata(finalSE_combo)$drug_combinations)==3)
print(all(sapply(metadata(finalSE_combo)$drug_combinations, 
    function(x) length(x$condition$Concentration_2) == (length(x$rows)-1))))
print(all(sapply(metadata(finalSE_combo)$drug_combinations, '[[', 'type') == 'fixed'))

write.table(df_merged_data, '../testdata/synthdata_combo_2dose_nonoise_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_combo, '../testdata/finalSE_combo_2dose_nonoise.RDS', compress = FALSE)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the test set with combo (two single dose)
#   co-treatment drug is also as single agent as DrugName
df_layout <- merge(CellLines[2:4,], Drugs[c(2:4,26),], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(Drugs[26*c(1,1,1),], Concentration = c(0, .2, 1))
colnames(df_2) <- paste0(colnames(df_2), '_2')
df_layout_2 <- merge(df_layout, df_2, by = NULL)

df_merged_data <- generate_response_data(df_layout_2, 0)
df_merged_data <- df_merged_data[!(df_merged_data$Gnumber %in% c('vehicle', Drugs$Gnumber[26]) & 
        df_merged_data$Gnumber_2 == Drugs$Gnumber[26]),]

finalSE_combo2 <- process_data_to_SE(df_merged_data)

# test accurarcy of the processing and fitting 
dt_test <- test_accuracy(finalSE_combo2[rowData(finalSE_combo2)$Concentration_2 == 0, ], e_inf, ec50, hill_coef)
print(dt_test)
# test: 
print(apply(abs(dt_test) < c(1e-3, 2e-3, 0.02, 0.015, 1e-4), 1, all))

# compare to other way of processing the data
DT1 <- convert_se_assay_to_dt(finalSE_combo, 'Metrics')
DT2 <- convert_se_assay_to_dt(finalSE_combo2,'Metrics')
DT1$rId <- gsub('_\\d\\d?$', '', DT1$rId)
DT2$rId <- gsub('_\\d\\d?$', '', DT2$rId)
# merge the two results
delta <- merge(DT1[ , c('rId', 'cId', 'normalization_type', 'x_0', 'x_max')], 
        DT2[ , c('rId', 'cId', 'normalization_type', 'x_0', 'x_max')], by = c('rId', 'cId', 'normalization_type'))
# test:
print(all(abs(quantile((delta$x_0.x - delta$x_0.y)[!grepl('vehicle', delta$rId)])) < .0005))

write.table(df_merged_data, '../testdata/synthdata_combo_2dose_nonoise2_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_combo2, '../testdata/finalSE_combo_2dose_nonoise2.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the 3rd test set with combo (two single dose)
#   co-treatment drug does NOT have single agent response
df_layout <- merge(CellLines[2:4,], Drugs[2:4,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(Drugs[c(26,26,26),], Concentration = c(0, .2, 1))
colnames(df_2) <- paste0(colnames(df_2), '_2')
df_layout_2 <- merge(df_layout, df_2, by = NULL)
df_layout_2 = df_layout_2[!(df_layout_2$Concentration == 0 & df_layout_2$Concentration_2 > 0), ]

df_merged_data <- generate_response_data(df_layout_2, 0)
finalSE_combo3 <- process_data_to_SE(df_merged_data)

dt_test <- test_accuracy(finalSE_combo3[rowData(finalSE_combo3)$Concentration_2 == 0, ], e_inf, ec50, hill_coef)
print(apply(abs(dt_test) < c(1e-3, 2e-3, 0.02, 0.015, 1e-4), 1, all))

# compare to the complete data
DT1 = convert_se_assay_to_dt(finalSE_combo,'Metrics')
DT3 = convert_se_assay_to_dt(finalSE_combo3,'Metrics')
DT1$rId = gsub('_\\d\\d?$', '', DT1$rId)
DT3$rId = gsub('_\\d\\d?$', '', DT3$rId)
delta = merge(DT1[ , c('rId', 'cId', 'normalization_type', 'x_0', 'x_inf', 'r2')], 
        DT3[ , c('rId', 'cId', 'normalization_type', 'x_0', 'x_inf', 'r2')], by = c('rId', 'cId', 'normalization_type'))

# checking the x_0 value was properly fitted for most cases (there are a few failures but it is ok)
print(sum(abs(delta$x_0.x - delta$x_0.y)>.004)<4)
print(sum(abs(delta$x_0.x - delta$x_0.y)>.025)<3)
# checking the x_inf value was properly fitted
sort(abs(delta$x_inf.x - delta$x_inf.y))
print(sum(abs(delta$x_inf.x - delta$x_inf.y)>.00003)<5)
print(sum(abs(delta$x_inf.x - delta$x_inf.y)>.008)<3)

write.table(df_merged_data, '../testdata/synthdata_combo_2dose_nonoise3_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_combo3, '../testdata/finalSE_combo_2dose_nonoise3.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the test set with combo (unique dose; many drug)
df_layout <- merge(CellLines[2:4,], Drugs[-1,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(Drugs[c(1,1),], Concentration = c(0, 2))
colnames(df_2) <- paste0(colnames(df_2), '_2')
df_layout_2 <- merge(df_layout, df_2, by = NULL)

df_merged_data <- generate_response_data(df_layout_2)
finalSE_combo <- process_data_to_SE(df_merged_data)

# test accuracy of the processing and fitting for the single agent
dt_test <- test_accuracy(finalSE_combo[rowData(finalSE_combo)$Concentration_2 == 0, ], e_inf, ec50, hill_coef)
print(apply(abs(dt_test) < c(0.5, 0.2, 2.5, 1.2, 0.3), 1, all))
# test the effect of the combination treatment
dt_test <- test_accuracy(finalSE_combo[rowData(finalSE_combo)$Concentration_2 > 0, ], e_inf, ec50, hill_coef)
print(apply(dt_test[c('delta_einf', '1_r2'),] < c(-.1, .01), 1, function(x) sum(x)==2))
# test the assignment of drug_combinations
print(length(metadata(finalSE_combo)$drug_combinations)==149)
print(all(sapply(metadata(finalSE_combo)$drug_combinations, 
    function(x) length(x$condition$Concentration_2) == (length(x$rows)-1))))
print(all(sapply(metadata(finalSE_combo)$drug_combinations, '[[', 'type') == 'fixed'))

write.table(df_merged_data, '../testdata/synthdata_combo_1dose_many_drugs_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_combo, '../testdata/finalSE_combo_1dose_many_drugs.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data with combo matrix (small, no noise)
df_layout <- merge(CellLines[7:8,], Drugs[c(4:6),], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout, Concentrations = 10**(seq(-3,.5,.5)))

df_2 <- merge(CellLines[CellLines$clid %in% df_layout$clid,], Drugs[c(21,26),], by = NULL)
df_2 <- add_data_replicates(df_2)
df_2 <- add_concentration(df_2, Concentrations = 10**(seq(-3,.5,.5)))
colnames(df_2)[colnames(df_2) %in% c(colnames(Drugs),'Concentration')] <- 
    paste0(colnames(df_2)[colnames(df_2) %in% c(colnames(Drugs),'Concentration')], '_2')

df_layout_2 <- merge(df_layout, df_2)

df_merged_data <- generate_response_data(df_layout_2, 0)
finalSE_matrix <- process_data_to_SE(df_merged_data)

# add and test calculation for combo matrix
# TODO when the functions are cleaned up


# test accuracy of the processing and fitting for the single agent
dt_test <- test_accuracy(finalSE_matrix[rowData(finalSE_matrix)$Concentration_2 == 0, ], e_inf, ec50, hill_coef)
print(apply(abs(dt_test) < c(1e-3, 6e-3, 0.12, 0.015, 1e-4), 1, all))
# test proper processing of the combo metadata
print(all(table(rowData(finalSE_matrix)[,c('Gnumber','Gnumber_2')])[, Drugs$Gnumber[c(21,26)]]==8))
print(all(table(rowData(finalSE_matrix)[rowData(finalSE_matrix)$DrugName_2 != 'vehicle','Concentration_2'])==6))

dt = convert_se_assay_to_dt(finalSE_matrix, 'Averaged')
print(all(dim(table(dt[dt$DrugName_2 != 'vehicle',c('Concentration', 'Concentration_2')]))==8))
print(all(table(dt[dt$DrugName_2 != 'vehicle',c('Concentration', 'Concentration_2')])==36))
# test the assignment of drug_combinations
print(length(metadata(finalSE_matrix)$drug_combinations)==6)
print(all(sapply(metadata(finalSE_matrix)$drug_combinations, 
    function(x) length(x$condition$Concentration_2) == (length(x$rows)-1))))
print(all(sapply(metadata(finalSE_matrix)$drug_combinations, '[[', 'type') == 'matrix'))

write.table(df_merged_data, '../testdata/synthdata_combo_matrix_small_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_matrix, '../testdata/finalSE_combo_matrix_small.RDS', compress = FALSE)


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data with combo matrix (mid-size)
df_layout <- merge(CellLines[seq(1,30,4),], Drugs[c(1,2,11,12,16,17),], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- merge(CellLines[CellLines$clid %in% df_layout$clid,], Drugs[c(21,26,31),], by = NULL)
df_2 <- add_data_replicates(df_2)
df_2 <- add_concentration(df_2)
colnames(df_2)[colnames(df_2) %in% c(colnames(Drugs),'Concentration')] <- 
    paste0(colnames(df_2)[colnames(df_2) %in% c(colnames(Drugs),'Concentration')], '_2')

df_layout_2 <- merge(df_layout, df_2)

df_merged_data <- generate_response_data(df_layout_2)
finalSE_matrix <- process_data_to_SE(df_merged_data)

# add and test calculation for combo matrix
# TODO when the functions are cleaned up

# test proper processing of the combo metadata
print(all(table(rowData(finalSE_matrix)[,c('Gnumber','Gnumber_2')])[, Drugs$Gnumber[c(21,26)]]==9))
print(all(table(rowData(finalSE_matrix)[rowData(finalSE_matrix)$DrugName_2 != 'vehicle','Concentration_2'])==18))

dt = convert_se_assay_to_dt(finalSE_matrix, 'Averaged')
print(all(dim(table(dt[dt$DrugName_2 != 'vehicle',c('Concentration', 'Concentration_2')]))==9))
print(all(table(dt[dt$DrugName_2 != 'vehicle',c('Concentration', 'Concentration_2')])==144))
# test the assignment of drug_combinations
print(length(metadata(finalSE_matrix)$drug_combinations)==18)
print(all(sapply(metadata(finalSE_matrix)$drug_combinations, 
    function(x) length(x$condition$Concentration_2) == (length(x$rows)-1))))
print(all(sapply(metadata(finalSE_matrix)$drug_combinations, '[[', 'type') == 'matrix'))

write.table(df_merged_data, '../testdata/synthdata_combo_matrix_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_matrix, '../testdata/finalSE_combo_matrix.RDS', compress = FALSE)



#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data with triple combo  (no noise)
df_layout <- merge(CellLines[7:8,], Drugs[c(4:6),], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout, Concentrations = 10**(seq(-3,.5,.5)))

df_2 <- merge(CellLines[CellLines$clid %in% df_layout$clid,], Drugs[c(21,26),], by = NULL)
df_2 <- add_data_replicates(df_2)
df_2 <- add_concentration(df_2, Concentrations = c(0, 10**(seq(-3,.5,.5))))
colnames(df_2)[colnames(df_2) %in% c(colnames(Drugs),'Concentration')] <- 
    paste0(colnames(df_2)[colnames(df_2) %in% c(colnames(Drugs),'Concentration')], '_2')

df_layout_2 <- merge(df_layout, df_2)

df_3 <- merge(CellLines[CellLines$clid %in% df_layout$clid,], Drugs[10,], by = NULL)
df_3 <- add_data_replicates(df_3)
df_3 <- add_concentration(df_3, Concentrations = c(0, .1, 1))
colnames(df_3)[colnames(df_3) %in% c(colnames(Drugs),'Concentration')] <- 
    paste0(colnames(df_3)[colnames(df_3) %in% c(colnames(Drugs),'Concentration')], '_3')

df_layout_3 <- merge(merge(df_layout, df_2), df_3)

df_merged_data <- generate_response_data(df_layout_3, 0)
finalSE_matrix <- process_data_to_SE(df_merged_data)

# add and test calculation for combo matrix
# TODO when the functions are cleaned up



# test accuracy of the processing and fitting for the single agent
dt_test <- test_accuracy(finalSE_matrix[rowData(finalSE_matrix)$Concentration_2 == 0 & 
    rowData(finalSE_matrix)$Concentration_3 == 0, ], e_inf, ec50, hill_coef)
print(apply(abs(dt_test) < c(1e-3, 6e-3, 0.12, 0.015, 1e-4), 1, all))
# test proper processing of the combo metadata
print(all(table(rowData(finalSE_matrix)[,c('Gnumber','Gnumber_2')])[, Drugs$Gnumber[c(21,26)]]==24))
print(all(table(rowData(finalSE_matrix)[rowData(finalSE_matrix)$DrugName_2 != 'vehicle','Concentration_2'])==18))
print(all(table(rowData(finalSE_matrix)[,paste0('Concentration',c('_2', '_3'))]) == c(3,array(6,8))))

dt = convert_se_assay_to_dt(finalSE_matrix, 'Averaged')
print(all(dim(table(dt[dt$DrugName_2 != 'vehicle',paste0('Concentration',c('', '_2', '_3'))]))==c(8,8,3)))
# test the assignment of drug_combinations
print(length(metadata(finalSE_matrix)$drug_combinations)==18)
print(all(sapply(metadata(finalSE_matrix)$drug_combinations, 
    function(x) length(x$condition$Concentration_2) == (length(x$rows)-1))))
print(all(sapply(metadata(finalSE_matrix)$drug_combinations, '[[', 'type') == 'matrix'))

write.table(df_merged_data, '../testdata/synthdata_combo_triple_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_matrix, '../testdata/finalSE_combo_triple.RDS', compress = FALSE)



#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data with combo co-dilution (small)
df_layout <- merge(CellLines[1:2,], Drugs[1:4,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(Drugs[1,,drop=F], df_layout[,'Concentration',drop=F])
colnames(df_2) <- paste0(colnames(df_2), '_2')
df_layout_2 <- cbind(df_layout, df_2)
df_layout_2 = df_layout_2[df_layout_2$DrugName != df_layout_2$DrugName_2,]
df_layout_2[df_layout_2$Concentration_2 > 0 , c('Concentration', 'Concentration_2')] <- 
  df_layout_2[df_layout_2$Concentration_2 > 0 , c('Concentration', 'Concentration_2')] / 2

df_merged_data <- generate_response_data(df_layout_2, 0)
finalSE_codilution <- process_data_to_SE(df_merged_data)

# add and test calculation for combo matrix
# TODO when the functions are cleaned up

write.table(df_merged_data, '../testdata/synthdata_combo_codilution_small_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_codilution, '../testdata/finalSE_combo_codilution_small.RDS', compress = FALSE)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# generate the data for the test set with combo (co-dilution)
df_layout <- merge(CellLines[seq(1,15,2),], Drugs[1:12,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(Drugs[c(1,1),], df_layout[,'Concentration',drop=F])
colnames(df_2) <- paste0(colnames(df_2), '_2')
df_layout_2 <- cbind(df_layout, df_2)
df_layout_2[df_layout_2$Concentration_2 > 0 , c('Concentration', 'Concentration_2')] <- 
  df_layout_2[df_layout_2$Concentration_2 > 0 , c('Concentration', 'Concentration_2')] / 2

df_merged_data <- generate_response_data(df_layout_2)
finalSE_codilution <- process_data_to_SE(df_merged_data)

# add and test calculation for combo matrix
# TODO when the functions are cleaned up

write.table(df_merged_data, '../testdata/synthdata_combo_codilution_rawdata.tsv', quote = F, row.names = F, sep = '\t')
saveRDS(finalSE_codilution, '../testdata/finalSE_combo_codilution.RDS', compress = FALSE)

