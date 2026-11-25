# prepareCodilutionData

Create data.table with input co-dilution data containing noise for
testing purposes

## Usage

``` r
prepareCodilutionData(
  cell_lines,
  drugs,
  drugsIdx2 = 1,
  conc = 10^(seq(-3, 1, 0.5)),
  noise = 0.1
)
```

## Arguments

- cell_lines:

  data.table with cell line info

- drugs:

  data.table with drug info

- drugsIdx2:

  numeric vector of ids for secondary drug (in `drugs` data.table)

- conc:

  vector of doses

- noise:

  number indicating level of noise

## Value

data.table with input data for testing

## Examples

``` r
prepareCodilutionData(create_synthetic_cell_lines()[seq_len(2), ],
                      create_synthetic_drugs()[seq_len(4), ])
#>      Barcode    clid CellLineName   Tissue ReferenceDivisionTime Gnumber
#>       <char>  <char>       <char>   <char>                 <num>  <char>
#>   1: plate_1 CL00010  cellline_AA tissue_x                    22 vehicle
#>   2: plate_1 CL00011  cellline_BA tissue_x                    26 vehicle
#>   3: plate_1 CL00010  cellline_AA tissue_x                    22 vehicle
#>   4: plate_1 CL00011  cellline_BA tissue_x                    26 vehicle
#>   5: plate_1 CL00010  cellline_AA tissue_x                    22 vehicle
#>  ---                                                                    
#> 194: plate_3 CL00011  cellline_BA tissue_x                    26  G00002
#> 195: plate_3 CL00010  cellline_AA tissue_x                    22  G00003
#> 196: plate_3 CL00011  cellline_BA tissue_x                    26  G00003
#> 197: plate_3 CL00010  cellline_AA tissue_x                    22  G00004
#> 198: plate_3 CL00011  cellline_BA tissue_x                    26  G00004
#>      DrugName drug_moa Concentration Gnumber_2 DrugName_2 drug_moa_2
#>        <char>   <char>         <num>    <char>     <char>     <char>
#>   1:  vehicle  vehicle             0   vehicle    vehicle    vehicle
#>   2:  vehicle  vehicle             0   vehicle    vehicle    vehicle
#>   3:  vehicle  vehicle             0   vehicle    vehicle    vehicle
#>   4:  vehicle  vehicle             0   vehicle    vehicle    vehicle
#>   5:  vehicle  vehicle             0   vehicle    vehicle    vehicle
#>  ---                                                                
#> 194: drug_002    moa_A             5    G00001   drug_001      moa_A
#> 195: drug_003    moa_A             5    G00001   drug_001      moa_A
#> 196: drug_003    moa_A             5    G00001   drug_001      moa_A
#> 197: drug_004    moa_A             5    G00001   drug_001      moa_A
#> 198: drug_004    moa_A             5    G00001   drug_001      moa_A
#>      Concentration_2 ReadoutValue BackgroundValue Duration
#>                <num>        <num>           <num>    <num>
#>   1:               0    101.90000               0       72
#>   2:               0     99.50000               0       72
#>   3:               0    101.60000               0       72
#>   4:               0    104.30000               0       72
#>   5:               0    104.40000               0       72
#>  ---                                                      
#> 194:               5     58.68946               0       72
#> 195:               5     22.30671               0       72
#> 196:               5     48.60656               0       72
#> 197:               5     38.94302               0       72
#> 198:               5     45.96414               0       72
```
