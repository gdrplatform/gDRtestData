# Copy Number Variation (CNV)

Gene-level copy number estimates from SNP microarray or WES.

## Format

Matrix with cell lines as rows and genes as columns

## Loading

`data.table::fread(system.file("depmap_data/OmicsCNGene.csv.gz", package = "gDRtestData"))`

## Description

- Rows: Cell line identifiers

- Columns: NCBI gene IDs

- Values: Log2 ratio relative to diploid reference

- Typical range: -2 (deletion) to +3 (amplification)

- 0 = diploid (2 copies)

## Source

[DepMap Portal - Data](https://depmap.org/portal/data_page/?tab=allData)

## Details

Downloaded May 26, 2026 from DepMap Portal (version 24Q4).

Citation: 24Q4 DepMap Release, including CRISPR Screens, PRISM Drug
Screens, Copy Number, Mutation, Expression, and Fusions DepMap, Broad
(2024). DepMap 24Q4 Public. Figshare+. Dataset.
<https://doi.org/10.25452/figshare.plus.27993248.v1>
