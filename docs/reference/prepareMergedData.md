# prepareMergedData

Create data.table with input data containing noise for testing purposes

## Usage

``` r
prepareMergedData(cell_lines, drugs, noise = 0.1)
```

## Arguments

- cell_lines:

  data.table with cell line info

- drugs:

  data.table with drug info

- noise:

  number indicating level of noise

## Value

data.table with input data for testing

## Examples

``` r
prepareMergedData(create_synthetic_cell_lines(), create_synthetic_drugs())
#>         Barcode    clid CellLineName   Tissue ReferenceDivisionTime Gnumber
#>          <char>  <char>       <char>   <char>                 <num>  <char>
#>      1: plate_1 CL00010  cellline_AA tissue_x                    22 vehicle
#>      2: plate_1 CL00011  cellline_BA tissue_x                    26 vehicle
#>      3: plate_1 CL00012  cellline_CA tissue_x                    30 vehicle
#>      4: plate_1 CL00013  cellline_DA tissue_x                    34 vehicle
#>      5: plate_1 CL00014  cellline_EA tissue_x                    38 vehicle
#>     ---                                                                    
#> 712796: plate_3 CL00095  cellline_KO tissue_z                    62  G00240
#> 712797: plate_3 CL00096  cellline_LO tissue_z                    66  G00240
#> 712798: plate_3 CL00097  cellline_MO tissue_z                    70  G00240
#> 712799: plate_3 CL00098  cellline_NO tissue_z                    74  G00240
#> 712800: plate_3 CL00099  cellline_OO tissue_z                    78  G00240
#>         DrugName drug_moa Concentration ReadoutValue BackgroundValue Duration
#>           <char>   <char>         <num>        <num>           <num>    <num>
#>      1:  vehicle  vehicle             0        100.7               0       72
#>      2:  vehicle  vehicle             0        103.6               0       72
#>      3:  vehicle  vehicle             0        104.6               0       72
#>      4:  vehicle  vehicle             0         98.7               0       72
#>      5:  vehicle  vehicle             0        102.4               0       72
#>     ---                                                                      
#> 712796: drug_240    moa_X            10         20.6               0       72
#> 712797: drug_240    moa_X            10          8.3               0       72
#> 712798: drug_240    moa_X            10         48.8               0       72
#> 712799: drug_240    moa_X            10         49.6               0       72
#> 712800: drug_240    moa_X            10         42.0               0       72
```
