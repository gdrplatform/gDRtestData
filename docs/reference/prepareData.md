# prepareData

Create data.table with input data for testing purposes

## Usage

``` r
prepareData(cell_lines, drugs, conc = 10^(seq(-3, 1, 0.5)))
```

## Arguments

- cell_lines:

  data.table with cell line info

- drugs:

  data.table with drug info

- conc:

  vector of doses

## Value

data.table with input data for testing

## Examples

``` r
prepareData(create_synthetic_cell_lines(), create_synthetic_drugs())
#>         Barcode    clid CellLineName   Tissue ReferenceDivisionTime Gnumber
#>          <char>  <char>       <char>   <char>                 <num>  <char>
#>      1: plate_1 CL00010  cellline_AA tissue_x                    22  G00001
#>      2: plate_1 CL00011  cellline_BA tissue_x                    26  G00001
#>      3: plate_1 CL00012  cellline_CA tissue_x                    30  G00001
#>      4: plate_1 CL00013  cellline_DA tissue_x                    34  G00001
#>      5: plate_1 CL00014  cellline_EA tissue_x                    38  G00001
#>     ---                                                                    
#> 712796: plate_3 CL00095  cellline_KO tissue_z                    62  G00240
#> 712797: plate_3 CL00096  cellline_LO tissue_z                    66  G00240
#> 712798: plate_3 CL00097  cellline_MO tissue_z                    70  G00240
#> 712799: plate_3 CL00098  cellline_NO tissue_z                    74  G00240
#> 712800: plate_3 CL00099  cellline_OO tissue_z                    78  G00240
#>         DrugName drug_moa Concentration
#>           <char>   <char>         <num>
#>      1: drug_001    moa_A             0
#>      2: drug_001    moa_A             0
#>      3: drug_001    moa_A             0
#>      4: drug_001    moa_A             0
#>      5: drug_001    moa_A             0
#>     ---                                
#> 712796: drug_240    moa_X            10
#> 712797: drug_240    moa_X            10
#> 712798: drug_240    moa_X            10
#> 712799: drug_240    moa_X            10
#> 712800: drug_240    moa_X            10
```
