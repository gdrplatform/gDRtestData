# ============================================================================
# Prepare DepMap 24Q4 Data for gDRtestData Package
# ============================================================================
#
# This script documents and reproduces the data preparation workflow.
#
# Source: DepMap Public 24Q4
# Downloaded: May 26, 2026
# URL: https://depmap.org/portal/data_page/?tab=allData
#
# Citation: DepMap, Broad (2024). DepMap 24Q4 Public. Figshare+. Dataset.
#           https://doi.org/10.25452/figshare.plus.27993248.v1
#
# ============================================================================

# ============================================================================
# SETUP
# ============================================================================

# Required packages
required_packages <- c("data.table", "withr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
}

source_dir <- "PATH_TO_DOWLOAD_DIR" # Update this path
final_dir <- "inst/depmap_data/"

# Verify directories exist
if (!dir.exists(source_dir)) {
  stop(sprintf("Source directory does not exist: %s", source_dir))
}

if (!dir.exists(final_dir)) {
  dir.create(final_dir, recursive = TRUE, showWarnings = TRUE)
  message(sprintf("Created output directory: %s", final_dir))
}

# ============================================================================
# 1. MODELS (Cell Line Metadata)
# ============================================================================
# Source file: Model.csv
file_name <- "Model.csv"

message("Processing Models metadata...")

depmap_models_raw <- data.table::fread(file.path(source_dir, file_name))

ls_col <- c(
  "ModelID",
  "CCLEName",
  "TissueOrigin",
  "DepmapModelType",
  "OncotreeLineage",
  "OncotreePrimaryDisease",
  "OncotreeSubtype",
  "OncotreeCode",
  "PrimaryOrMetastasis",
  "Age",
  "AgeCategory",
  "Sex",
  "PatientRace"
)

depmap_models <- depmap_models_raw[, .SD, .SDcols = ls_col]

write.csv(depmap_models,
  file = file.path(final_dir, file_name),
  row.names = FALSE
)

# ============================================================================
# 2. SOMATIC MUTATIONS - HOTSPOT
# ============================================================================
# Source file: OmicsSomaticMutationsMatrixHotspot.csv
file_name <- "OmicsSomaticMutationsMatrixHotspot.csv"

message("Processing Somatic Mutations (Hotspot)...")

mutations_hotspot_raw <- data.table::fread(file.path(source_dir, file_name))

id_col <- intersect(c("V1", "ModelID"), names(mutations_hotspot_raw))
ls_col_hotspot <- names(mutations_hotspot_raw)[names(mutations_hotspot_raw) != id_col]
ls_col_hotspot <- withr::with_seed(42, sample(ls_col_hotspot, 150))
depmap_mutations_hotspot <- mutations_hotspot_raw[, .SD, .SDcols = c(id_col, ls_col_hotspot)]

write.csv(depmap_mutations_hotspot,
  file = file.path(final_dir, file_name),
  row.names = FALSE
)

# ============================================================================
# 3. SOMATIC MUTATIONS - DAMAGING
# ============================================================================
# Source file: OmicsSomaticMutationsMatrixDamaging.csv
file_name <- "OmicsSomaticMutationsMatrixDamaging.csv"

message("Processing Somatic Mutations (Damaging)...")

mutations_damaging_raw <- data.table::fread(file.path(source_dir, file_name))

id_col <- intersect(c("V1", "ModelID"), names(mutations_damaging_raw))
ls_col <- intersect(names(mutations_damaging_raw), ls_col_hotspot)
depmap_mutations_damaging <- mutations_damaging_raw[, .SD, .SDcols = c(id_col, ls_col)]

write.csv(depmap_mutations_damaging,
  file = file.path(final_dir, file_name),
  row.names = FALSE
)

# ============================================================================
# 4. CRISPR GENE EFFECT
# ============================================================================
# Source file: CRISPRGeneEffect.csv
file_name <- "CRISPRGeneEffect.csv"

message("Processing CRISPR Gene Effect...")

crispr_raw <- data.table::fread(file.path(source_dir, file_name))

id_col <- intersect(c("V1", "ModelID"), names(crispr_raw))
ls_col <- intersect(names(crispr_raw), ls_col_hotspot)
depmap_crispr_gene_effect <- crispr_raw[, .SD, .SDcols = c(id_col, ls_col)]
depmap_crispr_gene_effect[, (ls_col) := lapply(.SD, round, 4), .SDcols = ls_col]

write.csv(depmap_crispr_gene_effect,
  file = file.path(final_dir, file_name),
  row.names = FALSE
)

# ============================================================================
# 5. GENE EXPRESSION (Protein-Coding Genes, TPM, Log-transformed)
# ============================================================================
# Source file: OmicsExpressionProteinCodingGenesTPMLogp1.csv
file_name <- "OmicsExpressionProteinCodingGenesTPMLogp1.csv"

message("Processing Gene Expression...")

expression_raw <- data.table::fread(file.path(source_dir, file_name))

id_col <- intersect(c("V1", "ModelID"), names(expression_raw))
ls_col <- intersect(names(expression_raw), ls_col_hotspot)
depmap_expression <- expression_raw[, .SD, .SDcols = c(id_col, ls_col)]
depmap_expression[, (ls_col) := lapply(.SD, round, 4), .SDcols = ls_col]

write.csv(depmap_expression,
  file = file.path(final_dir, file_name),
  row.names = FALSE
)

# ============================================================================
# 6. COPY NUMBER VARIATION (CNV)
# ============================================================================
# Source file: OmicsCNGene.csv
file_name <- "OmicsCNGene.csv"

message("Processing Copy Number Variation...")

cn_raw <- data.table::fread(file.path(source_dir, file_name))

id_col <- intersect(c("V1", "ModelID"), names(cn_raw))
ls_col <- intersect(names(cn_raw), ls_col_hotspot)
depmap_copy_number <- cn_raw[, .SD, .SDcols = c(id_col, ls_col)]
depmap_copy_number[, (ls_col) := lapply(.SD, round, 4), .SDcols = ls_col]

write.csv(depmap_copy_number,
  file = file.path(final_dir, file_name),
  row.names = FALSE
)

# ============================================================================
# 7. VALIDATION
# ============================================================================

message("\n=== Data Consistency Check ===")
# Ensure all datasets have consistent cell line IDs
common_models <- Reduce(
  intersect,
  list(
    depmap_models[[1]],
    depmap_crispr_gene_effect[[1]],
    depmap_expression[[1]],
    depmap_mutations_hotspot[[1]],
    depmap_mutations_damaging[[1]],
    depmap_copy_number[[1]]
  )
)
cat(sprintf("Common models across all datasets (%d models)", NROW(common_models)))


# Ensure all datasets have consistent genes
common_genes <- Reduce(
  intersect,
  list(
    names(depmap_crispr_gene_effect)[-1],
    names(depmap_expression)[-1],
    names(depmap_mutations_hotspot)[-1],
    names(depmap_mutations_damaging)[-1],
    names(depmap_copy_number)[-1]
  )
)
cat(sprintf("Common genes across all datasets (%d genes)", NROW(common_genes)))

# ============================================================================
# 8. DATASET INFORMATION
# ============================================================================

dataset_info <- list(
  version = "24Q4",
  download_date = "2026-05-26",
  source_url = "https://depmap.org/portal/data_page/?tab=allData",
  citation = "DepMap, Broad (2024). DepMap 24Q4 Public. Figshare+. https://doi.org/10.25452/figshare.plus.27993248.v1",
  datasets = list(
    models = list(
      rows = NROW(depmap_models),
      cols = NCOL(depmap_models)
    ),
    crispr_gene_effect = list(
      rows = NROW(depmap_crispr_gene_effect),
      cols = NCOL(depmap_crispr_gene_effect)
    ),
    expression = list(
      rows = NROW(depmap_expression),
      cols = NCOL(depmap_expression)
    ),
    mutations_hotspot = list(
      rows = NROW(depmap_mutations_hotspot),
      cols = NCOL(depmap_mutations_hotspot)
    ),
    mutations_damaging = list(
      rows = NROW(depmap_mutations_damaging),
      cols = NCOL(depmap_mutations_damaging)
    ),
    copy_number = list(
      rows = NROW(depmap_copy_number),
      cols = NCOL(depmap_copy_number)
    )
  )
)
