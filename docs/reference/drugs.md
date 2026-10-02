# Drugs

Drugs

## Value

data.table

## Examples

``` r
path <- system.file("annotation_data", "drugs.csv", package = "gDRtestData")
data.table::fread(file = path)
#>      Gnumber    DrugName drug_moa
#>       <char>      <char>   <char>
#>   1: vehicle     vehicle  vehicle
#>   2:  G00002    drug_002    moa_A
#>   3:  G00003    drug_003    moa_A
#>   4:  G00004    drug_004    moa_A
#>   5:  G00005    drug_005    moa_A
#>  ---                             
#> 758:  G00844 PF-05175157    ACACA
#> 759:  G00845  BMS-303141     ACLY
#> 760:  G00846   LY3295668    AURKA
#> 761:  G00847   FASN-IN-4     FASN
#> 762:  G00848   Sotorasib     KRAS
```
