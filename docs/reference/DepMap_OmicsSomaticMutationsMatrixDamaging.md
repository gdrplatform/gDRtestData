# Somatic Mutations (Damaging)

Binary matrix of damaging mutations (frame-shift, stop-gain,
splice-site).

## Format

Matrix with cell lines as rows and genes as columns

## Loading

`data.table::fread(system.file("depmap_data/OmicsSomaticMutationsMatrixDamaging.csv.gz", package = "gDRtestData"))`

## Description

- Rows: Cell line identifiers

- Columns: NCBI gene IDs

- Values: Binary (0 = no mutation, 1 = damaging mutation)

- High-confidence loss-of-function mutations

## Source

[DepMap Portal - Data](https://depmap.org/portal/data_page/?tab=allData)

## Details

Downloaded May 26, 2026 from DepMap Portal (version 24Q4).

Citation: 24Q4 DepMap Release, including CRISPR Screens, PRISM Drug
Screens, Copy Number, Mutation, Expression, and Fusions DepMap, Broad
(2024). DepMap 24Q4 Public. Figshare+. Dataset.
<https://doi.org/10.25452/figshare.plus.27993248.v1>
