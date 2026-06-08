# CRISPR Gene Effect Scores

Genome-wide CRISPR/Cas9 knockout dependency scores.

## Format

Matrix with cell lines as rows and genes as columns

## Loading

`data.table::fread(system.file("depmap_data/CRISPRGeneEffect.csv.gz", package = "gDRtestData"))`

## Description

- Rows: Cell line identifiers

- Columns: NCBI gene IDs (Entrez format)

- Values: Dependency scores (-1 to +1); lower = more essential

- NA indicates insufficient screen coverage

## Source

[DepMap Portal - Data](https://depmap.org/portal/data_page/?tab=allData)

## Details

Downloaded May 26, 2026 from DepMap Portal (version 24Q4).

Citation: 24Q4 DepMap Release, including CRISPR Screens, PRISM Drug
Screens, Copy Number, Mutation, Expression, and Fusions DepMap, Broad
(2024). DepMap 24Q4 Public. Figshare+. Dataset.
<https://doi.org/10.25452/figshare.plus.27993248.v1>
