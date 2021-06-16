# datasets for testing purposes only:
#   Cases that could raise issues with data processing
#
#  1. No barcode
#  2. random masking
#  3. one cell line masked
#  4. one drug masked


library(SummarizedExperiment)
library(BumpyMatrix)
library(gDRutils)
library(gDRwrapper)
library(gDRcore)
library(reshape2)

cell_lines <- create_synthetic_cell_lines()
drugs <- create_synthetic_drugs()
e_inf <- generate_e_inf(drugs, cell_lines)
ec50 <- generate_ec50(drugs, cell_lines)
hill_coef <- generate_hill_coef(drugs, cell_lines)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#### test "barcode"
# generate the data without barcode
df_layout <- merge(cell_lines[2:5,], drugs[2:6,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout, 0)

# no barcode
df_merged_data <- df_merged_data
df_merged_data$Barcode <- NULL

finalSE_1_no_noise <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# test accurarcy of the processing and fitting (no noise => low tolerance)
dt_test <- test_accuracy(finalSE_1_no_noise, e_inf, ec50, hill_coef)
# test:
print(apply(abs(dt_test) < c(1e-3, 2.2e-3, 0.05, 0.015, 1e-4), 1, all))



#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#### test missing data
# generate the data for the 1st test set: no noise
#   without masked
df_layout <- merge(cell_lines[2:5,], drugs[2:6,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout, 0)

finalSE_1_no_noise <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# test accurarcy of the processing and fitting (no noise => low tolerance)
dt_test <- test_accuracy(finalSE_1_no_noise, e_inf, ec50, hill_coef)
# test:
print(apply(abs(dt_test) < c(1e-3, 2.2e-3, 0.05, 0.015, 1e-4), 1, all))


# remove some samples by cell line/drug
df_merged_data_missing <- df_merged_data[df_merged_data$clid != "CL00011", ]

finalSE_1_no_noise <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data_missing)

# test accurarcy of the processing and fitting (no noise => low tolerance)
dt_test <- test_accuracy(finalSE_1_no_noise, e_inf, ec50, hill_coef)
# test:
print(apply(abs(dt_test) < c(1e-3, 5e-3, 0.05, 0.015, 1e-4), 1, all))



df_merged_data_missing <- df_merged_data[df_merged_data$clid != "CL00012" &
            df_merged_data$Gnumber != "G00004", ]
finalSE_1_no_noise <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data_missing)


df_merged_data_missing <- df_merged_data[df_merged_data$clid != "CL00011" |
            df_merged_data$Gnumber != "G00002", ]
finalSE_1_no_noise <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data_missing)



#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#### test "masked"
# generate the data for the 1st test set: no noise
#   with masked = F
df_layout <- merge(cell_lines[2:5,], drugs[2:6,], by = NULL)
df_layout <- add_data_replicates(df_layout)
df_layout <- add_concentration(df_layout)

df_merged_data <- generate_response_data(df_layout, 0)

df_merged_data$masked <- FALSE

finalSE_1_no_noise <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data)

# test accurarcy of the processing and fitting (no noise => low tolerance)
dt_test <- test_accuracy(finalSE_1_no_noise, e_inf, ec50, hill_coef)
# test:
print(apply(abs(dt_test) < c(1e-3, 2.2e-3, 0.05, 0.015, 1e-4), 1, all))


# mask some samples randomly
df_merged_data <- generate_response_data(df_layout, 0)
df_merged_data_masked <- df_merged_data
set.seed(2)
df_merged_data_masked$masked <- runif(nrow(df_merged_data))>.9

finalSE_1_no_noise <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data_masked)

# test accurarcy of the processing and fitting (no noise => low tolerance)
dt_test <- test_accuracy(finalSE_1_no_noise, e_inf, ec50, hill_coef)
# test:
print(apply(abs(dt_test) < c(1e-3, 5e-3, 0.05, 0.015, 1e-4), 1, all))


# mask some samples randomly (40%)
df_merged_data <- generate_response_data(df_layout, 0)
df_merged_data_masked <- df_merged_data
set.seed(2)
df_merged_data_masked$masked <- runif(nrow(df_merged_data)) > .6

finalSE_1_no_noise <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data_masked)

# test accurarcy of the processing and fitting (no noise => low tolerance)
dt_test <- test_accuracy(finalSE_1_no_noise, e_inf, ec50, hill_coef)
# test:
print(apply(abs(dt_test) < c(1e-3, 5e-3, 0.06, 0.015, 1e-4), 1, all))



# mask some samples for one cell line/drug
df_merged_data <- generate_response_data(df_layout, 0)
df_merged_data_masked <- df_merged_data
df_merged_data_masked$masked <- FALSE
df_merged_data_masked$masked[df_merged_data$clid == "CL00011" & df_merged_data$Gnumber == "G00004"] <- TRUE


finalSE_1_no_noise <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data_masked)

print(nrow(assay(finalSE_1_no_noise, "Averaged")[rowData(finalSE_1_no_noise)$Gnumber == "G00004",
    colData(finalSE_1_no_noise)$clid == "CL00011"]) == 1)

print(is.na(assay(finalSE_1_no_noise, "Averaged")[rowData(finalSE_1_no_noise)$Gnumber == "G00004",
    colData(finalSE_1_no_noise)$clid == "CL00011"][[1]][1,"RelativeViability"]))


# mask all samples of a cell line
df_merged_data <- generate_response_data(df_layout, 0)
df_merged_data_masked = df_merged_data
df_merged_data_masked$masked = F
df_merged_data_masked$masked[df_merged_data$clid == "CL00011"] = TRUE

finalSE_1_no_noise <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data_masked)


se <- gDR::create_SE2(df_merged_data_masked)
normSE <- gDR::normalize_SE2(se)
avgSE <- gDR::average_SE2(normSE)
metricsSE <- gDR::fit_SE2(avgSE)
finalSE <- gDR::add_codrug_group_SE(metricsSE)


print(all(sapply(assay(finalSE_1_no_noise, "Averaged")[,
    colData(finalSE_1_no_noise)$clid == "CL00011"], nrow)==1))

print(all(sapply(assay(finalSE_1_no_noise, "Averaged")[,
    colData(finalSE_1_no_noise)$clid == "CL00011"],
    function(x) is.na(x[1,"RelativeViability"]))))

print(all(sapply(assay(finalSE_1_no_noise, "Metrics")[,
    colData(finalSE_1_no_noise)$clid == "CL00011"],
    function(x) is.na(x["RV", "x_mean"]))))


# mask all samples of a drug
df_merged_data <- generate_response_data(df_layout, 0)
df_merged_data_masked = df_merged_data
df_merged_data_masked$masked = F
df_merged_data_masked$masked[df_merged_data$Gnumber == "G00004"] = TRUE

finalSE_1_no_noise <- gDRcore::runDrugResponseProcessingPipeline(df_merged_data_masked)


print(all(sapply(assay(finalSE_1_no_noise, "Averaged")[
    rowData(finalSE_1_no_noise)$Gnumber == "G00004", ], nrow)==1))

print(all(sapply(assay(finalSE_1_no_noise, "Averaged")[
    rowData(finalSE_1_no_noise)$Gnumber == "G00004", ],
    function(x) is.na(x[1,"RelativeViability"]))))

print(all(sapply(assay(finalSE_1_no_noise, "Metrics")[
    rowData(finalSE_1_no_noise)$Gnumber == "G00004", ],
    function(x) is.na(x["RV", "x_mean"]))))

