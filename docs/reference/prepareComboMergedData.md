# prepareComboMergedData

Create data.table with input combination data containing noise for
testing purposes

## Usage

``` r
prepareComboMergedData(
  cell_lines,
  drugs,
  drugsIdx1 = 2:4,
  drugsIdx2 = c(26, 26, 26),
  concentration = c(0, 0.2, 1),
  noise = 0.1,
  modifyDf2 = FALSE
)
```

## Arguments

- cell_lines:

  data.table with cell line info

- drugs:

  data.table with drug info

- drugsIdx1:

  numeric vector of ids for primary drug

- drugsIdx2:

  numeric vector of ids for secondary drug

- concentration:

  numeric vector of doses

- noise:

  number indicating level of noise

- modifyDf2:

  Boolean indicating if the table should me modified to keep reverse
  single agent data

## Value

data.table with input data for testing

## Examples

``` r
prepareComboMergedData(create_synthetic_cell_lines(), create_synthetic_drugs())
#>        Barcode    clid CellLineName   Tissue ReferenceDivisionTime Gnumber
#>         <char>  <char>       <char>   <char>                 <num>  <char>
#>     1: plate_1 CL00010  cellline_AA tissue_x                    22 vehicle
#>     2: plate_1 CL00011  cellline_BA tissue_x                    26 vehicle
#>     3: plate_1 CL00012  cellline_CA tissue_x                    30 vehicle
#>     4: plate_1 CL00013  cellline_DA tissue_x                    34 vehicle
#>     5: plate_1 CL00014  cellline_EA tissue_x                    38 vehicle
#>    ---                                                                    
#> 26726: plate_3 CL00095  cellline_KO tissue_z                    62  G00004
#> 26727: plate_3 CL00096  cellline_LO tissue_z                    66  G00004
#> 26728: plate_3 CL00097  cellline_MO tissue_z                    70  G00004
#> 26729: plate_3 CL00098  cellline_NO tissue_z                    74  G00004
#> 26730: plate_3 CL00099  cellline_OO tissue_z                    78  G00004
#>        DrugName drug_moa Concentration Gnumber_2 DrugName_2 drug_moa_2
#>          <char>   <char>         <num>    <char>     <char>     <char>
#>     1:  vehicle  vehicle             0   vehicle    vehicle    vehicle
#>     2:  vehicle  vehicle             0   vehicle    vehicle    vehicle
#>     3:  vehicle  vehicle             0   vehicle    vehicle    vehicle
#>     4:  vehicle  vehicle             0   vehicle    vehicle    vehicle
#>     5:  vehicle  vehicle             0   vehicle    vehicle    vehicle
#>    ---                                                                
#> 26726: drug_004    moa_A            10    G00026   drug_026      moa_E
#> 26727: drug_004    moa_A            10    G00026   drug_026      moa_E
#> 26728: drug_004    moa_A            10    G00026   drug_026      moa_E
#> 26729: drug_004    moa_A            10    G00026   drug_026      moa_E
#> 26730: drug_004    moa_A            10    G00026   drug_026      moa_E
#>        Concentration_2 ReadoutValue BackgroundValue Duration
#>                  <num>        <num>           <num>    <num>
#>     1:               0    99.900000               0       72
#>     2:               0    99.800000               0       72
#>     3:               0   103.100000               0       72
#>     4:               0    97.800000               0       72
#>     5:               0   104.700000               0       72
#>    ---                                                      
#> 26726:               1     3.301784               0       72
#> 26727:               1     8.842755               0       72
#> 26728:               1     2.313932               0       72
#> 26729:               1    10.121007               0       72
#> 26730:               1    25.243907               0       72
```
