# Drugs

Drugs

## Value

data.table

## Examples

``` r
path <- system.file("annotation_data", "drugs.csv", package = "gDRtestData")
data.table::fread(file = path)
#>     Gnumber     DrugName                drug_moa
#>      <char>       <char>                  <char>
#>  1: vehicle      vehicle                 vehicle
#>  2:  G00002     drug_002                   moa_A
#>  3:  G00003     drug_003                   moa_A
#>  4:  G00004     drug_004                   moa_A
#>  5:  G00005     drug_005                   moa_A
#>  6:  G00006     drug_006                   moa_A
#>  7:  G00007     drug_007                   moa_A
#>  8:  G00008     drug_008                   moa_A
#>  9:  G00009     drug_009                   moa_A
#> 10:  G00010     drug_010                   moa_A
#> 11:  G00011     drug_011                   moa_B
#> 12:  G00021     drug_021                   moa_D
#> 13:  G00026     drug_026                   moa_E
#> 14:  G00100   Everolimus                    MTOR
#> 15:  G00101   Inavolisib                  PIK3CA
#> 16:  G00102 Giredestrant                    ESR1
#> 17:  G00103 Belvarafenib          BRAF|RAF1|ARAF
#> 18:  G00104  Vemurafenib                    BRAF
#> 19:  G00105  Cobimetinib           MAP2K1|MAP2K2
#> 20:  G00106     GDC-8025 TEAD1|TEAD2|TEAD3|TEAD4
#>     Gnumber     DrugName                drug_moa
#>      <char>       <char>                  <char>
```
