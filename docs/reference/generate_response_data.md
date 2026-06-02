# Generate response data

Generate response data

## Usage

``` r
generate_response_data(df_layout, noise_level = 0.1)
```

## Arguments

- df_layout:

  data.table that should contains the cell line, drug, concentration,
  and replicate columns along with the annotations that needs to be
  propagated

- noise_level:

  numeric scalar with the level of noise added to the data

## Value

data.table with response data

## Examples

``` r
cell_lines <- create_synthetic_cell_lines()
drugs <- create_synthetic_drugs()
df_layout <- prepareData(cell_lines[seq_len(2), ], drugs[seq_len(4), ])
generate_response_data(df_layout)
#>      Barcode    clid CellLineName   Tissue ReferenceDivisionTime Gnumber
#>       <char>  <char>       <char>   <char>                 <num>  <char>
#>   1: plate_1 CL00010  cellline_AA tissue_x                    22 vehicle
#>   2: plate_1 CL00011  cellline_BA tissue_x                    26 vehicle
#>   3: plate_1 CL00010  cellline_AA tissue_x                    22 vehicle
#>   4: plate_1 CL00011  cellline_BA tissue_x                    26 vehicle
#>   5: plate_1 CL00010  cellline_AA tissue_x                    22 vehicle
#>  ---                                                                    
#> 260: plate_3 CL00011  cellline_BA tissue_x                    26  G00002
#> 261: plate_3 CL00010  cellline_AA tissue_x                    22  G00003
#> 262: plate_3 CL00011  cellline_BA tissue_x                    26  G00003
#> 263: plate_3 CL00010  cellline_AA tissue_x                    22  G00004
#> 264: plate_3 CL00011  cellline_BA tissue_x                    26  G00004
#>      DrugName drug_moa Concentration ReadoutValue BackgroundValue Duration
#>        <char>   <char>         <num>        <num>           <num>    <num>
#>   1:  vehicle  vehicle             0        101.2               0       72
#>   2:  vehicle  vehicle             0         96.2               0       72
#>   3:  vehicle  vehicle             0         99.6               0       72
#>   4:  vehicle  vehicle             0         97.6               0       72
#>   5:  vehicle  vehicle             0        101.6               0       72
#>  ---                                                                      
#> 260: drug_002    moa_A            10         75.9               0       72
#> 261: drug_003    moa_A            10         81.3               0       72
#> 262: drug_003    moa_A            10         73.7               0       72
#> 263: drug_004    moa_A            10         68.9               0       72
#> 264: drug_004    moa_A            10         71.1               0       72

```
