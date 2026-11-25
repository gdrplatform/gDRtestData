# Cell lines

Cell lines

## Value

data.table

## Examples

``` r
path <- system.file("annotation_data", "cell_lines.csv", package = "gDRtestData")
data.table::fread(file = path)
#>        clid CellLineName Tissue ReferenceDivisionTime parental_identifier
#>      <char>       <char> <char>                 <int>              <char>
#>  1: CL00016  cellline_GB breast                    46         cellline_GB
#>  2: CL00017  cellline_HB breast                    50         cellline_HB
#>  3: CL00011  cellline_BA breast                    26         cellline_BA
#>  4: CL00012  cellline_CA breast                    30         cellline_CA
#>  5: CL00013  cellline_DA breast                    34         cellline_DA
#>  6: CL00014  cellline_EA breast                    38         cellline_EA
#>  7: CL00015  cellline_FA breast                    42         cellline_FA
#>  8: CL00018  cellline_IB breast                    54         cellline_IB
#>  9: CL00019  cellline_JB breast                    58         cellline_JB
#> 10: CL00020  cellline_KB breast                    62         cellline_KB
#>     subtype
#>      <lgcl>
#>  1:      NA
#>  2:      NA
#>  3:      NA
#>  4:      NA
#>  5:      NA
#>  6:      NA
#>  7:      NA
#>  8:      NA
#>  9:      NA
#> 10:      NA
```
