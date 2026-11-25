# Add data replicates

Add data replicates

## Usage

``` r
add_data_replicates(df_layout)
```

## Arguments

- df_layout:

  data.table that should contains the cell line, drug, concentration,
  and replicate columns along with the annotations that needs to be
  propagated

## Value

data.table with replicates

## Examples

``` r
cell_lines <- create_synthetic_cell_lines()
add_data_replicates(cell_lines)
#>      Barcode    clid CellLineName   Tissue ReferenceDivisionTime
#>       <char>  <char>       <char>   <char>                 <num>
#>   1: plate_1 CL00010  cellline_AA tissue_x                    22
#>   2: plate_1 CL00011  cellline_BA tissue_x                    26
#>   3: plate_1 CL00012  cellline_CA tissue_x                    30
#>   4: plate_1 CL00013  cellline_DA tissue_x                    34
#>   5: plate_1 CL00014  cellline_EA tissue_x                    38
#>  ---                                                            
#> 266: plate_3 CL00095  cellline_KO tissue_z                    62
#> 267: plate_3 CL00096  cellline_LO tissue_z                    66
#> 268: plate_3 CL00097  cellline_MO tissue_z                    70
#> 269: plate_3 CL00098  cellline_NO tissue_z                    74
#> 270: plate_3 CL00099  cellline_OO tissue_z                    78
```
