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
library(BumpyMatrix)
library(gDRutils)
library(gDRwrapper)
library(gDRcore)
library(reshape2)


# Helper functions
save_tsv <- function(object, filename) {
  write.table(object, system.file("testdata", filename, package = "gDRtestData"), quote = FALSE, row.names = FALSE, sep = "\t")
}

save_rds <- function(object, filename) {
  saveRDS(object, system.file("testdata", filename, package = "gDRtestData"), compress = FALSE)
}

cell_lines <- create_synthetic_cell_lines()
drugs <- create_synthetic_drugs()
e_inf <- generate_e_inf(drugs, cell_lines)
ec50 <- generate_ec50(drugs, cell_lines)
hill_coef <- generate_hill_coef(drugs, cell_lines)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 1st test set: no noise
#   only for testing purpuses not displayed as example
df_layout <- merge(cell_lines[2:11,], drugs[2:11,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout, 0)

finalMAE_1_no_noise <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# test accurarcy of the processing and fitting (no noise => low tolerance)
dt_test <- test_accuracy(finalMAE_1_no_noise[[1]], e_inf, ec50, hill_coef)
print(dt_test)
# test:
print(apply(abs(dt_test) < c(1e-3, 2.2e-3, 0.04, 0.015, 1e-4), 1, all))


save_tsv(df_merged_data, "synthdata_small_no_noise_rawdata.tsv")
save_rds(finalMAE_1_no_noise, "finalMAE_small_no_noise.RDS")

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 1st test set with noise
df_layout <- merge(cell_lines[2:11,], drugs[2:11,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout)
finalMAE_1 <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# test accurarcy of the processing and fitting (noise => medium tolerance)
dt_test <- test_accuracy(finalMAE_1[[1]], e_inf, ec50, hill_coef)
print(dt_test)
# test:
print(apply(abs(dt_test) < c(0.5, 0.1, 1.5, 1.2, 0.05), 1, all))

save_tsv(df_merged_data, "synthdata_small_rawdata.tsv")
save_rds(finalMAE_1_no_noise, "finalMAE_small.RDS")

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 1st test set with ligand as reference
df_layout <- merge(cell_lines[2:6,], drugs[2:5,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout, 0)
df_merged_data$Ligand <- 0.1
df_merged_data2 <- df_merged_data[df_merged_data$Gnumber %in% c("vehicle", "G00002", "G00003"),]
df_merged_data2$Ligand <- 0
df_merged_data2$ReadoutValue <- 105 - pmax(0, pmin(104, (105-df_merged_data2$ReadoutValue) ^ 1.1))
df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0("CL000", 11:12)] <-
    0.8 * df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0("CL000", 11:12)]
df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0("CL000", 13:14)] <-
    0.5 * df_merged_data2$ReadoutValue[df_merged_data2$clid %in% paste0("CL000", 13:14)]
df_merged_data2$ReadoutValue <- round(df_merged_data2$ReadoutValue,1)
df_merged_data2$Barcode <- paste0(df_merged_data2$Barcode, "1")
df_merged_data <- rbind(df_merged_data, df_merged_data2)

finalMAE_1_Ligand <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data, override_untrt_controls = c(Ligand = 0.1))

# test accurarcy of the processing and fitting for Ligand = 0.1 (noise => medium tolerance)
dt_test <- test_accuracy(finalMAE_1_Ligand[[1]][rowData(finalMAE_1_Ligand[[1]])$Ligand > 0,], e_inf, ec50, hill_coef)
print(dt_test)
# test:
print(apply(abs(dt_test) < c(1e-3, 3e-3, 0.031, 0.015, 1e-4), 1, all))
# test fit quality for Ligand = 0 and that delta(e_inf) < 0
dt_test <- test_accuracy(finalMAE_1_Ligand[[1]][rowData(finalMAE_1_Ligand[[1]])$Ligand == 0, ], e_inf, ec50, hill_coef)
print(dt_test)
# test:
print(apply(dt_test[c("delta_einf", "1_r2"),] < c(-0.15, 1e-4), 1, all))

save_tsv(df_merged_data, "finalSE_wLigand_rawdata.tsv")
save_rds(finalMAE_1_Ligand, "finalMAE_wLigand.RDS")


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 2nd (medium size) test set with single agent
df_layout <- merge(cell_lines[1:15,], drugs[1:40,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout)
finalMAE_2 <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# test accurarcy of the processing and fitting
dt_test <- test_accuracy(finalMAE_2[[1]], e_inf, ec50, hill_coef)
print(dt_test)
# test:
print(apply(abs(dt_test) < c(0.5, 0.2, 2.5, 1.2, 0.3), 1, all))

save_tsv(df_merged_data, "synthdata_medium_rawdata.tsv")
save_rds(finalMAE_2, "finalMAE_medium.RDS")


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 3rd (many lines) test set with single agent
df_layout <- merge(cell_lines, drugs[1:40,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout)
finalMAE_3 <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# test accurarcy of the processing and fitting
dt_test <- test_accuracy(finalMAE_3[[1]], e_inf, ec50, hill_coef)
print(dt_test)
# test:
print(apply(abs(dt_test) < c(0.5, 0.2, 2.5, 1.2, 0.3), 1, all))

save_tsv(df_merged_data, "synthdata_many_lines_rawdata.tsv")
save_rds(finalMAE_1_Ligand, "finalMAE_many_lines.RDS")


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the test set with single agent (many drugs)
df_layout <- merge(cell_lines[1:10,], drugs, by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout)
finalMAE_4 <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# test accurarcy of the processing and fitting
dt_test <- test_accuracy(finalMAE_4[[1]], e_inf, ec50, hill_coef)
print(dt_test)
# test:
print(apply(abs(dt_test) < c(0.5, 0.2, 2.5, 1.2, 0.3), 1, all))

save_tsv(df_merged_data, "synthdata_many_drugs_rawdata.tsv")
save_rds(finalMAE_4, "finalMAE_many_drugs.RDS")


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the test set with combo (two single dose)
#   co-treatment drug is only as DrugName_2
df_layout <- merge(cell_lines[2:4,], drugs[2:4,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(drugs[c(26,26,26),], Concentration = c(0, .2, 1))
colnames(df_2) <- paste0(colnames(df_2), "_2")
df_layout_2 <- merge(df_layout, df_2, by = NULL)

df_merged_data <- generate_response_data(df_layout_2, 0)
finalMAE_combo <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# test the single agent response
dt_test <- test_accuracy(finalMAE_combo[[1]][rowData(finalMAE_combo[[1]])$Concentration_2 == 0, ], e_inf, ec50, hill_coef)
print(dt_test)
# test:
print(apply(abs(dt_test) < c(1e-3, 2e-3, 0.02, 0.015, 1e-4), 1, all))
# test the assignment of drug_combinations
print(length(metadata(finalMAE_combo[[1]])$drug_combinations)==3)
print(all(sapply(metadata(finalMAE_combo[[1]])$drug_combinations,
    function(x) length(x$condition$Concentration_2) == (length(x$rows)-1))))
print(all(sapply(metadata(finalSE_combo[[1]])$drug_combinations, "[[", "type") == "fixed"))

save_tsv(df_merged_data, "synthdata_combo_2dose_nonoise_rawdata.tsv")
save_rds(finalMAE_combo, "finalMAE_combo_2dose_nonoise.RDS")

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the test set with combo (two single dose)
#   co-treatment drug is also as single agent as DrugName
df_layout <- merge(cell_lines[2:4,], drugs[c(2:4,26),], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(drugs[26*c(1,1,1),], Concentration = c(0, .2, 1))
colnames(df_2) <- paste0(colnames(df_2), "_2")
df_layout_2 <- merge(df_layout, df_2, by = NULL)

df_merged_data <- generate_response_data(df_layout_2, 0)
df_merged_data <- df_merged_data[!(df_merged_data$Gnumber %in% c("vehicle", drugs$Gnumber[26]) &
        df_merged_data$Gnumber_2 == drugs$Gnumber[26]),]

finalMAE_combo2 <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# test accurarcy of the processing and fitting
dt_test <- test_accuracy(finalMAE_combo2[[1]][rowData(finalMAE_combo2[[1]])$Concentration_2 == 0, ], e_inf, ec50, hill_coef)
print(dt_test)
# test:
print(apply(abs(dt_test) < c(1e-3, 2e-3, 0.02, 0.015, 1e-4), 1, all))

# compare to other way of processing the data
DT1 <- convert_se_assay_to_dt(finalMAE_combo[[1]], "Metrics")
DT2 <- convert_se_assay_to_dt(finalMAE_combo2[[1]],"Metrics")
DT1$rId <- gsub("_\\d\\d?$", "", DT1$rId)
DT2$rId <- gsub("_\\d\\d?$", "", DT2$rId)
# merge the two results
delta <- merge(DT1[ , c("rId", "cId", "normalization_type", "x_0", "x_max")],
        DT2[ , c("rId", "cId", "normalization_type", "x_0", "x_max")], by = c("rId", "cId", "normalization_type"))
# test:
print(all(abs(quantile((delta$x_0.x - delta$x_0.y)[!grepl("vehicle", delta$rId)])) < .0005))

save_tsv(df_merged_data, "synthdata_combo_2dose_nonoise2_rawdata.tsv")
save_rds(finalMAE_combo2, "finalMAE_combo_2dose_nonoise2.RDS")

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 3rd test set with combo (two single dose)
#   co-treatment drug does NOT have single agent response
df_layout <- merge(cell_lines[2:4,], drugs[2:4,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(drugs[c(26,26,26),], Concentration = c(0, .2, 1))
colnames(df_2) <- paste0(colnames(df_2), "_2")
df_layout_2 <- merge(df_layout, df_2, by = NULL)
df_layout_2 = df_layout_2[!(df_layout_2$Concentration == 0 & df_layout_2$Concentration_2 > 0), ]

df_merged_data <- generate_response_data(df_layout_2, 0)
finalMAE_combo3 <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

dt_test <- test_accuracy(finalMAE_combo3[[1]][rowData(finalMAE_combo3[[1]])$Concentration_2 == 0, ], e_inf, ec50, hill_coef)
print(apply(abs(dt_test) < c(1e-3, 2e-3, 0.02, 0.015, 1e-4), 1, all))

# compare to the complete data
DT1 = convert_se_assay_to_dt(finalMAE_combo[[1]],"Metrics")
DT3 = convert_se_assay_to_dt(finalMAE_combo3[[1]],"Metrics")
DT1$rId = gsub("_\\d\\d?$", "", DT1$rId)
DT3$rId = gsub("_\\d\\d?$", "", DT3$rId)
delta = merge(DT1[ , c("rId", "cId", "normalization_type", "x_0", "x_inf", "r2")],
        DT3[ , c("rId", "cId", "normalization_type", "x_0", "x_inf", "r2")], by = c("rId", "cId", "normalization_type"))

# checking the x_0 value was properly fitted for most cases (there are a few failures but it is ok)
print(sum(abs(delta$x_0.x - delta$x_0.y)>.004)<4)
print(sum(abs(delta$x_0.x - delta$x_0.y)>.025)<3)
# checking the x_inf value was properly fitted
sort(abs(delta$x_inf.x - delta$x_inf.y))
print(sum(abs(delta$x_inf.x - delta$x_inf.y)>.00003)<5)
print(sum(abs(delta$x_inf.x - delta$x_inf.y)>.008)<3)

save_tsv(df_merged_data, "synthdata_combo_2dose_nonoise3_rawdata.tsv")
save_rds(finalMAE_combo3, "finalMAE_combo_2dose_nonoise3.RDS")

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the test set with combo (unique dose; many drug)
df_layout <- merge(cell_lines[2:4,], drugs[-1,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(drugs[c(1,1),], Concentration = c(0, 2))
colnames(df_2) <- paste0(colnames(df_2), "_2")
df_layout_2 <- merge(df_layout, df_2, by = NULL)

df_merged_data <- generate_response_data(df_layout_2)
finalMAE_combo <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# test accuracy of the processing and fitting for the single agent
dt_test <- test_accuracy(finalMAE_combo[[1]][rowData(finalMAE_combo[[1]])$Concentration_2 == 0, ], e_inf, ec50, hill_coef)
print(apply(abs(dt_test) < c(0.5, 0.2, 2.5, 1.2, 0.3), 1, all))
# test the effect of the combination treatment
dt_test <- test_accuracy(finalMAE_combo[[1]][rowData(finalMAE_combo[[1]])$Concentration_2 > 0, ], e_inf, ec50, hill_coef)
print(apply(dt_test[c("delta_einf", "1_r2"),] < c(-.1, .01), 1, function(x) sum(x)==2))
# test the assignment of drug_combinations
print(length(metadata(finalMAE_combo)$drug_combinations)==149)
print(all(sapply(metadata(finalMAE_combo[[1]])$drug_combinations,
    function(x) length(x$condition$Concentration_2) == (length(x$rows)-1))))
print(all(sapply(metadata(finalMAE_combo[[1]])$drug_combinations, "[[", "type") == "fixed"))

save_tsv(df_merged_data, "synthdata_combo_1dose_many_drugs_rawdata.tsv")
save_rds(finalMAE_combo, "finalMAE_combo_1dose_many_drugs.RDS")

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data with combo matrix (small, no noise)
df_layout <- merge(cell_lines[7:8,], drugs[c(4:6),], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout, concentrations = 10^ (seq(-3,.5,.5)))

df_2 <- merge(cell_lines[cell_lines$clid %in% df_layout$clid,], drugs[c(21,26),], by = NULL)
df_2 <- add_data_replicates(df_2)
df_2 <- add_concentration(df_2, concentrations = 10^ (seq(-3,.5,.5)))
colnames(df_2)[colnames(df_2) %in% c(colnames(drugs),"Concentration")] <-
    paste0(colnames(df_2)[colnames(df_2) %in% c(colnames(drugs),"Concentration")], "_2")

df_layout_2 <- merge(df_layout, df_2)

df_merged_data <- generate_response_data(df_layout_2, 0)
finalMAE_matrix <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# add and test calculation for combo matrix
# TODO when the functions are cleaned up


# test accuracy of the processing and fitting for the single agent
dt_test <- test_accuracy(finalMAE_matrix[[1]][rowData(finalMAE_matrix[[1]])$Concentration_2 == 0, ], e_inf, ec50, hill_coef)
print(apply(abs(dt_test) < c(1e-3, 6e-3, 0.12, 0.015, 1e-4), 1, all))
# test proper processing of the combo metadata
print(all(table(rowData(finalMAE_matrix[[1]])[,c("Gnumber","Gnumber_2")])[, drugs$Gnumber[c(21,26)]]==8))
print(all(table(rowData(finalMAE_matrix[[1]])[rowData(finalMAE_matrix[[1]])$DrugName_2 != "vehicle","Concentration_2"])==6))

dt = convert_se_assay_to_dt(finalMAE_matrix[[1]], "Averaged")
print(all(dim(table(dt[dt$DrugName_2 != "vehicle",c("Concentration", "Concentration_2")]))==8))
print(all(table(dt[dt$DrugName_2 != "vehicle",c("Concentration", "Concentration_2")])==36))
# test the assignment of drug_combinations
print(length(metadata(finalMAE_matrix[[1]])$drug_combinations)==6)
print(all(sapply(metadata(finalMAE_matrix[[1]])$drug_combinations,
    function(x) length(x$condition$Concentration_2) == (length(x$rows)-1))))
print(all(sapply(metadata(finalMAE_matrix[[1]])$drug_combinations, "[[", "type") == "matrix"))

save_tsv(df_merged_data, "synthdata_combo_matrix_small_rawdata.tsv")
save_rds(finalMAE_matrix, "finalMAE_combo_matrix_small.RDS")


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data with combo matrix (mid-size)
df_layout <- merge(cell_lines[seq(1,30,4),], drugs[c(1,2,11,12,16,17),], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- merge(cell_lines[cell_lines$clid %in% df_layout$clid,], drugs[c(21,26,31),], by = NULL)
df_2 <- add_data_replicates(df_2)
df_2 <- add_concentration(df_2)
colnames(df_2)[colnames(df_2) %in% c(colnames(drugs),"Concentration")] <-
    paste0(colnames(df_2)[colnames(df_2) %in% c(colnames(drugs),"Concentration")], "_2")

df_layout_2 <- merge(df_layout, df_2)

df_merged_data <- generate_response_data(df_layout_2)
finalMAE_matrix <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# add and test calculation for combo matrix
# TODO when the functions are cleaned up

# test proper processing of the combo metadata
print(all(table(rowData(finalMAE_matrix[[1]])[,c("Gnumber","Gnumber_2")])[, drugs$Gnumber[c(21,26)]]==9))
print(all(table(rowData(finalMAE_matrix[[1]])[rowData(finalMAE_matrix[[1]])$DrugName_2 != "vehicle","Concentration_2"])==18))

dt = convert_se_assay_to_dt(finalMAE_matrix[[1]], "Averaged")
print(all(dim(table(dt[dt$DrugName_2 != "vehicle",c("Concentration", "Concentration_2")]))==9))
print(all(table(dt[dt$DrugName_2 != "vehicle",c("Concentration", "Concentration_2")])==144))
# test the assignment of drug_combinations
print(length(metadata(finalMAE_matrix[[1]])$drug_combinations)==18)
print(all(sapply(metadata(finalMAE_matrix[[1]])$drug_combinations,
    function(x) length(x$condition$Concentration_2) == (length(x$rows)-1))))
print(all(sapply(metadata(finalMAE_matrix[[1]])$drug_combinations, "[[", "type") == "matrix"))

save_tsv(df_merged_data, "synthdata_combo_matrix_rawdata.tsv")
save_rds(finalMAE_combo2, "finalMAE_combo_matrix.RDS")


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data with triple combo  (no noise)
df_layout <- merge(cell_lines[7:8,], drugs[c(4:6),], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout, concentrations = 10^ (seq(-3,.5,.5)))

df_2 <- merge(cell_lines[cell_lines$clid %in% df_layout$clid,], drugs[c(21,26),], by = NULL)
df_2 <- add_data_replicates(df_2)
df_2 <- add_concentration(df_2, concentrations = c(0, 10^ (seq(-3,.5,.5))))
colnames(df_2)[colnames(df_2) %in% c(colnames(drugs),"Concentration")] <-
    paste0(colnames(df_2)[colnames(df_2) %in% c(colnames(drugs),"Concentration")], "_2")

df_layout_2 <- merge(df_layout, df_2)

df_3 <- merge(cell_lines[cell_lines$clid %in% df_layout$clid,], drugs[10,], by = NULL)
df_3 <- add_data_replicates(df_3)
df_3 <- add_concentration(df_3, concentrations = c(0, .1, 1))
colnames(df_3)[colnames(df_3) %in% c(colnames(drugs),"Concentration")] <-
    paste0(colnames(df_3)[colnames(df_3) %in% c(colnames(drugs),"Concentration")], "_3")

df_layout_3 <- merge(merge(df_layout, df_2), df_3)

df_merged_data <- generate_response_data(df_layout_3, 0)
finalMAE_matrix <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# add and test calculation for combo matrix
# TODO when the functions are cleaned up



# test accuracy of the processing and fitting for the single agent
dt_test <- test_accuracy(finalMAE_matrix[[1]][rowData(finalMAE_matrix[[1]])$Concentration_2 == 0 &
    rowData(finalMAE_matrix[[1]])$Concentration_3 == 0, ], e_inf, ec50, hill_coef)
print(apply(abs(dt_test) < c(1e-3, 6e-3, 0.12, 0.015, 1e-4), 1, all))
# test proper processing of the combo metadata
print(all(table(rowData(finalMAE_matrix[[1]])[,c("Gnumber","Gnumber_2")])[, drugs$Gnumber[c(21,26)]]==24))
print(all(table(rowData(finalMAE_matrix[[1]])[rowData(finalMAE_matrix[[1]])$DrugName_2 != "vehicle","Concentration_2"])==18))
print(all(table(rowData(finalMAE_matrix[[1]])[,paste0("Concentration",c("_2", "_3"))]) == c(3,array(6,8))))

dt = convert_se_assay_to_dt(finalMAE_matrix[[1]], "Averaged")
print(all(dim(table(dt[dt$DrugName_2 != "vehicle",paste0("Concentration",c("", "_2", "_3"))]))==c(8,8,3)))
# test the assignment of drug_combinations
print(length(metadata(finalMAE_matrix[[1]])$drug_combinations)==18)
print(all(sapply(metadata(finalMAE_matrix[[1]])$drug_combinations,
    function(x) length(x$condition$Concentration_2) == (length(x$rows)-1))))
print(all(sapply(metadata(finalMAE_matrix[[1]])$drug_combinations, "[[", "type") == "matrix"))

save_tsv(df_merged_data, "synthdata_combo_triple_rawdata.tsv")
save_rds(finalMAE_matrix, "finalMAE_combo_triple.RDS")


#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data with combo co-dilution (small)
df_layout <- merge(cell_lines[1:2,], drugs[1:4,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(drugs[1,,drop= FALSE], df_layout[,"Concentration",drop= FALSE])
colnames(df_2) <- paste0(colnames(df_2), "_2")
df_layout_2 <- cbind(df_layout, df_2)
df_layout_2 = df_layout_2[df_layout_2$DrugName != df_layout_2$DrugName_2,]
df_layout_2[df_layout_2$Concentration_2 > 0 , c("Concentration", "Concentration_2")] <-
  df_layout_2[df_layout_2$Concentration_2 > 0 , c("Concentration", "Concentration_2")] / 2

df_merged_data <- generate_response_data(df_layout_2, 0)
finalMAE_codilution <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# add and test calculation for combo matrix
# TODO when the functions are cleaned up

save_tsv(df_merged_data, "synthdata_combo_codilution_small_rawdata.tsv")
save_rds(finalMAE_codilution, "finalMAE_combo_codilution_small.RDS")

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the test set with combo (co-dilution)
df_layout <- merge(cell_lines[seq(1,15,2),], drugs[1:12,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_2 <- cbind(drugs[c(1,1),], df_layout[,"Concentration",drop= FALSE])
colnames(df_2) <- paste0(colnames(df_2), "_2")
df_layout_2 <- cbind(df_layout, df_2)
df_layout_2[df_layout_2$Concentration_2 > 0 , c("Concentration", "Concentration_2")] <-
  df_layout_2[df_layout_2$Concentration_2 > 0 , c("Concentration", "Concentration_2")] / 2

df_merged_data <- generate_response_data(df_layout_2)
finalMAE_codilution <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# add and test calculation for combo matrix
# TODO when the functions are cleaned up

save_tsv(df_merged_data, "synthdata_combo_codilution_rawdata.tsv")
save_rds(finalMAE_codilution, "finalMAE_combo_codilution.RDS")
