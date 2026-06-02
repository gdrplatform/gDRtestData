# Add data with day 0

Add data with day 0

## Usage

``` r
add_day0_data(df_merged, noise_level = 0.05)
```

## Arguments

- df_merged:

  data.table with merged data

- noise_level:

  numeric scalar with the level of noise added to the data

## Value

data.table with day0 data

## Examples

``` r
cell_lines <- create_synthetic_cell_lines()
drugs <- create_synthetic_drugs()
df_merged <- prepareData(cell_lines[seq_len(2), ], drugs[seq_len(4), ])
df_merged$Duration <- 72
df_merged$ReadoutValue <- 0
add_day0_data(df_merged)
#>      Barcode    clid CellLineName   Tissue ReferenceDivisionTime Gnumber
#>       <char>  <char>       <char>   <char>                 <num>  <char>
#>   1: plate_1 CL00010  cellline_AA tissue_x                    22  G00001
#>   2: plate_1 CL00011  cellline_BA tissue_x                    26  G00001
#>   3: plate_1 CL00010  cellline_AA tissue_x                    22  G00002
#>   4: plate_1 CL00011  cellline_BA tissue_x                    26  G00002
#>   5: plate_1 CL00010  cellline_AA tissue_x                    22  G00003
#>  ---                                                                    
#> 284: plate_0 CL00011  cellline_BA tissue_x                    26  G00002
#> 285: plate_0 CL00010  cellline_AA tissue_x                    22  G00003
#> 286: plate_0 CL00011  cellline_BA tissue_x                    26  G00003
#> 287: plate_0 CL00010  cellline_AA tissue_x                    22  G00004
#> 288: plate_0 CL00011  cellline_BA tissue_x                    26  G00004
#>      DrugName drug_moa Concentration Duration ReadoutValue
#>        <char>   <char>         <num>    <num>        <num>
#>   1: drug_001    moa_A             0       72            0
#>   2: drug_001    moa_A             0       72            0
#>   3: drug_002    moa_A             0       72            0
#>   4: drug_002    moa_A             0       72            0
#>   5: drug_003    moa_A             0       72            0
#>  ---                                                      
#> 284: drug_002    moa_A             0        0            0
#> 285: drug_003    moa_A             0        0            0
#> 286: drug_003    moa_A             0        0            0
#> 287: drug_004    moa_A             0        0            0
#> 288: drug_004    moa_A             0        0            0
```
