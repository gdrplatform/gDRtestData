# generating three types of example/test data
#
#  The data have no day0 information.
#  dataset with * are to be imported as example for visualization
#
#  "small_no_noise"           10 drugs (3 differents drug_moa) by 10 lines (3 tissues); 
#                             single agent - no noise in the data
#  "small"                *    10 drugs (3 differents drug_moa) by 10 lines (3 tissues); single agent
#  "wLigand"              *    3 drugs by 4 lines (3 tissues); "Ligand = 0.1" as reference; single agent
#  "medium"               *    40 drugs (6 differents drug_moa) by 15 lines (3 tissues); single agent
#  "many_lines"           *    150 drugs (6 differents drug_moa) by 10 lines (3 tissues); single agent
#  "many_drugs"           *    150 drugs (6 differents drug_moa) by 10 lines (3 tissues); single agent
#  "combo_2dose_nonoise"  *    3 drugs x 2 co-treatment (1 drug at 2 doses) by 3 lines; 
#                              co-treatment drug is only as DrugName_2
#  "combo_2dose_nonoise2"      3 drugs x 2 co-treatment (1 drug at 2 doses) by 3 lines; 
#                              co-treatment drug is also as single agent as DrugName
#  "combo_2dose_nonoise3"      3 drugs x 2 co-treatment (1 drug at 2 doses) by 3 lines; 
#                              co-treatment drug does NOT have single agent response
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

set.seed(2)

cell_lines <- create_synthetic_cell_lines()
drugs <- create_synthetic_drugs()
e_inf <- generate_e_inf(drugs, cell_lines)
ec50 <- generate_ec50(drugs, cell_lines)
hill_coef <- generate_hill_coef(drugs, cell_lines)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 1st test set: no noise
#   only for testing purpuses not displayed as example
generateNoNoiseRawData(cell_lines, drugs, e_inf, ec50, hill_coef)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 1st test set with noise
generateNoiseRawData(cell_lines, drugs, e_inf, ec50, hill_coef)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 1st test set with ligand as reference
generateLigandData(cell_lines, drugs, e_inf, ec50, hill_coef)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 2nd (medium size) test set with single agent
generateMediumData(cell_lines, drugs, e_inf, ec50, hill_coef)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 3rd (many lines) test set with single agent
generateManyLinesData(cell_lines, drugs, e_inf, ec50, hill_coef)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the test set with single agent (many drugs)
generateManyDrugsData(cell_lines, drugs, e_inf, ec50, hill_coef)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the test set with combo (two single dose)
#   co-treatment drug is only as DrugName_2
generateComboNoNoiseData(cell_lines, drugs, e_inf, ec50, hill_coef)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the test set with combo (two single dose)
#   co-treatment drug is also as single agent as DrugName
generateComboNoNoiseData2(cell_lines, drugs, e_inf, ec50, hill_coef)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 3rd test set with combo (two single dose)
#   co-treatment drug does NOT have single agent response
generateComboNoNoiseData3(cell_lines, drugs, e_inf, ec50, hill_coef)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the test set with combo (unique dose; many drug)
generateComboManyDrugs(cell_lines, drugs, e_inf, ec50, hill_coef)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data with combo matrix (small, no noise)
generateComboMatrixSmall(cell_lines, drugs, e_inf, ec50, hill_coef)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data with combo matrix (mid-size)
generateComboMatrix(cell_lines, drugs, e_inf, ec50, hill_coef)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data with triple combo  (no noise)
generateTripleComboMatrix(cell_lines, drugs, e_inf, ec50, hill_coef)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data with combo co-dilution (small)
generateCodilutionSmall(cell_lines, drugs, e_inf, ec50, hill_coef)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the test set with combo (co-dilution)
generateCodilution(cell_lines, drugs, e_inf, ec50, hill_coef)

