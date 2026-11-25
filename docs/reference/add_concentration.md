# Add concentrations

Add concentrations

## Usage

``` r
add_concentration(df_layout, concentrations = 10^(seq(-3, 1, 0.5)))
```

## Arguments

- df_layout:

  data.table that should contains the cell line, drug, concentration,
  and replicate columns along with the annotations that needs to be
  propagated

- concentrations:

  vector of numeric concentrations that will be added to df_layout

## Value

data.table with concentrations

## Examples

``` r
cell_lines <- create_synthetic_cell_lines()
add_concentration(cell_lines)
#>         clid CellLineName   Tissue ReferenceDivisionTime Concentration
#>       <char>       <char>   <char>                 <num>         <num>
#>   1: CL00010  cellline_AA tissue_x                    22             0
#>   2: CL00011  cellline_BA tissue_x                    26             0
#>   3: CL00012  cellline_CA tissue_x                    30             0
#>   4: CL00013  cellline_DA tissue_x                    34             0
#>   5: CL00014  cellline_EA tissue_x                    38             0
#>  ---                                                                  
#> 986: CL00095  cellline_KO tissue_z                    62            10
#> 987: CL00096  cellline_LO tissue_z                    66            10
#> 988: CL00097  cellline_MO tissue_z                    70            10
#> 989: CL00098  cellline_NO tissue_z                    74            10
#> 990: CL00099  cellline_OO tissue_z                    78            10
```
