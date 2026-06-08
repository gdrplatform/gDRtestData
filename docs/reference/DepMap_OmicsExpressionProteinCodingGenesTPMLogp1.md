# Gene Expression (Log2 TPM)

RNA-seq based gene expression for protein-coding genes.

## Format

Matrix with cell lines as rows and genes as columns

## Loading

`data.table::fread( system.file("depmap_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv.gz", package = "gDRtestData"))`

## Description

- Rows: Cell line identifiers

- Columns: NCBI gene IDs (Entrez format)

- Values: Log2(TPM + 1) transformed expression

- Only protein-coding genes included

## Source

[DepMap Portal - Data](https://depmap.org/portal/data_page/?tab=allData)

## Details

Downloaded May 26, 2026 from DepMap Portal (version 24Q4).

Citation: 24Q4 DepMap Release, including CRISPR Screens, PRISM Drug
Screens, Copy Number, Mutation, Expression, and Fusions DepMap, Broad
(2024). DepMap 24Q4 Public. Figshare+. Dataset.
<https://doi.org/10.25452/figshare.plus.27993248.v1>
