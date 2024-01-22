# generating three types of example/test data
#
#  The data have no day0 information.
#  dataset with * are to be imported as example for visualization
#
#  "small_no_noise"            10 drugs (3 different drug_moa) by 10 lines (3 tissues); 
#                              single agent - no noise in the data
#  "small"                  *  10 drugs (3 different drug_moa) by 10 lines (3 tissues); single agent
#  "wLigand"                *  3 drugs by 4 lines (3 tissues); "Ligand = 0.1" as reference; single agent
#  "medium"                 *  40 drugs (6 different drug_moa) by 15 lines (3 tissues); single agent
#  "many_lines"             *  150 drugs (6 different drug_moa) by 10 lines (3 tissues); single agent
#  "many_drugs"             *  150 drugs (6 different drug_moa) by 10 lines (3 tissues); single agent
#  "combo_2dose_nonoise"    *  3 drugs x 2 co-treatment (1 drug at 2 doses) by 3 cell lines; 
#                              co-treatment drug occurs also as a primary drug
#  "combo_2dose_nonoise2"      3 drugs x 2 co-treatment (1 drug at 2 doses) by 3 cell lines; 
#                              co-treatment drug is also as single agent as DrugName
#  "combo_2dose_nonoise3"      3 drugs x 2 co-treatment (1 drug at 2 doses) by 3 cell lines; 
#                              co-treatment drug does NOT have single agent response
#  "combo_1dose_many_drugs" *  149 drugs x 1 drug (1 dose) by 3 lines;
#  "combo_matrix_small"        3 x 2 drugs (matrix) for 2 cell lines; no noise
#  "combo_matrix"           *  6 x 3 drugs (matrix) for 8 cell lines
#  "combo_triple"           *  2 x 3 x 2 drugs (triple matrix) for 4 cell lines
#  "combo_codilution_small"    4 x 1 drugs (co-dilution) for 2 cell lines; no noise
#  "combo_codilution"          12 x 1 drugs (co-dilution) for 8 cell lines; no noise

library(SummarizedExperiment)
library(BumpyMatrix)
library(gDRutils)
library(gDRcore)


set.seed(2)

cell_lines <- create_synthetic_cell_lines()
drugs <- create_synthetic_drugs()

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 1st test set: no noise
#   only for testing purpuses not displayed as example
set.seed(2)
generateNoNoiseRawData(cell_lines, drugs)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 1st test set with noise
set.seed(2)
generateNoiseRawData(cell_lines, drugs)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 1st test set with ligand as reference
set.seed(2)
generateLigandData(cell_lines, drugs)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 2nd (medium size) test set with single agent
set.seed(2)
generateMediumData(cell_lines, drugs)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the test set with combo (two single dose)
#   co-treatment drug is only as DrugName_2
set.seed(2)
generateComboNoNoiseData(cell_lines, drugs)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the test set with combo (two single dose)
#   co-treatment drug is also as single agent as DrugName
set.seed(2)
generateComboNoNoiseData2(cell_lines, drugs)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the 3rd test set with combo (two single dose)
#   co-treatment drug does NOT have single agent response
set.seed(2)
generateComboNoNoiseData3(cell_lines, drugs)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data with combo matrix (small, no noise)
set.seed(2)
generateComboMatrixSmall(cell_lines, drugs)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data with combo matrix (mid-size)
set.seed(2)
generateComboMatrix(cell_lines, drugs)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data with triple combo  (no noise)
set.seed(2)
generateTripleComboMatrix(cell_lines, drugs)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data with combo co-dilution (small)
set.seed(2)
generateCodilutionSmall(cell_lines, drugs)

#### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# generate the data for the test set with combo (co-dilution)
set.seed(2)
generateCodilution(cell_lines, drugs)

