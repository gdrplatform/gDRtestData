# Create data.table with synthetic cell lines

Create data.table with synthetic cell lines

## Usage

``` r
create_synthetic_cell_lines()
```

## Value

data.table with synthetic cell lines

## Examples

``` r
create_synthetic_cell_lines()
#>        clid CellLineName   Tissue ReferenceDivisionTime
#>      <char>       <char>   <char>                 <num>
#>  1: CL00010  cellline_AA tissue_x                    22
#>  2: CL00011  cellline_BA tissue_x                    26
#>  3: CL00012  cellline_CA tissue_x                    30
#>  4: CL00013  cellline_DA tissue_x                    34
#>  5: CL00014  cellline_EA tissue_x                    38
#>  6: CL00015  cellline_FA tissue_x                    42
#>  7: CL00016  cellline_GB tissue_y                    46
#>  8: CL00017  cellline_HB tissue_y                    50
#>  9: CL00018  cellline_IB tissue_y                    54
#> 10: CL00019  cellline_JB tissue_z                    58
#> 11: CL00020  cellline_KB tissue_z                    62
#> 12: CL00021  cellline_LB tissue_z                    66
#> 13: CL00022  cellline_MC tissue_z                    70
#> 14: CL00023  cellline_NC tissue_z                    74
#> 15: CL00024  cellline_OC tissue_z                    78
#> 16: CL00025  cellline_AC tissue_w                    22
#> 17: CL00026  cellline_BC tissue_w                    26
#> 18: CL00027  cellline_CC tissue_w                    30
#> 19: CL00028  cellline_DD tissue_w                    34
#> 20: CL00029  cellline_ED tissue_w                    38
#> 21: CL00030  cellline_FD tissue_w                    42
#> 22: CL00031  cellline_GD tissue_w                    46
#> 23: CL00032  cellline_HD tissue_w                    50
#> 24: CL00033  cellline_ID tissue_w                    54
#> 25: CL00034  cellline_JE tissue_w                    58
#> 26: CL00035  cellline_KE tissue_w                    62
#> 27: CL00036  cellline_LE tissue_w                    66
#> 28: CL00037  cellline_ME tissue_w                    70
#> 29: CL00038  cellline_NE tissue_w                    74
#> 30: CL00039  cellline_OE tissue_w                    78
#> 31: CL00040  cellline_AF tissue_w                    22
#> 32: CL00041  cellline_BF tissue_w                    26
#> 33: CL00042  cellline_CF tissue_w                    30
#> 34: CL00043  cellline_DF tissue_w                    34
#> 35: CL00044  cellline_EF tissue_w                    38
#> 36: CL00045  cellline_FF tissue_w                    42
#> 37: CL00046  cellline_GG tissue_w                    46
#> 38: CL00047  cellline_HG tissue_w                    50
#> 39: CL00048  cellline_IG tissue_w                    54
#> 40: CL00049  cellline_JG tissue_w                    58
#> 41: CL00050  cellline_KG tissue_v                    62
#> 42: CL00051  cellline_LG tissue_v                    66
#> 43: CL00052  cellline_MH tissue_v                    70
#> 44: CL00053  cellline_NH tissue_v                    74
#> 45: CL00054  cellline_OH tissue_v                    78
#> 46: CL00055  cellline_AH tissue_v                    22
#> 47: CL00056  cellline_BH tissue_v                    26
#> 48: CL00057  cellline_CH tissue_v                    30
#> 49: CL00058  cellline_DI tissue_v                    34
#> 50: CL00059  cellline_EI tissue_v                    38
#> 51: CL00060  cellline_FI tissue_x                    42
#> 52: CL00061  cellline_GI tissue_y                    46
#> 53: CL00062  cellline_HI tissue_y                    50
#> 54: CL00063  cellline_II tissue_y                    54
#> 55: CL00064  cellline_JJ tissue_z                    58
#> 56: CL00065  cellline_KJ tissue_z                    62
#> 57: CL00066  cellline_LJ tissue_z                    66
#> 58: CL00067  cellline_MJ tissue_z                    70
#> 59: CL00068  cellline_NJ tissue_z                    74
#> 60: CL00069  cellline_OJ tissue_z                    78
#> 61: CL00070  cellline_AK tissue_x                    22
#> 62: CL00071  cellline_BK tissue_x                    26
#> 63: CL00072  cellline_CK tissue_x                    30
#> 64: CL00073  cellline_DK tissue_x                    34
#> 65: CL00074  cellline_EK tissue_x                    38
#> 66: CL00075  cellline_FK tissue_x                    42
#> 67: CL00076  cellline_GL tissue_y                    46
#> 68: CL00077  cellline_HL tissue_y                    50
#> 69: CL00078  cellline_IL tissue_y                    54
#> 70: CL00079  cellline_JL tissue_z                    58
#> 71: CL00080  cellline_KL tissue_z                    62
#> 72: CL00081  cellline_LL tissue_z                    66
#> 73: CL00082  cellline_MM tissue_z                    70
#> 74: CL00083  cellline_NM tissue_z                    74
#> 75: CL00084  cellline_OM tissue_z                    78
#> 76: CL00085  cellline_AM tissue_x                    22
#> 77: CL00086  cellline_BM tissue_x                    26
#> 78: CL00087  cellline_CM tissue_x                    30
#> 79: CL00088  cellline_DN tissue_x                    34
#> 80: CL00089  cellline_EN tissue_x                    38
#> 81: CL00090  cellline_FN tissue_x                    42
#> 82: CL00091  cellline_GN tissue_y                    46
#> 83: CL00092  cellline_HN tissue_y                    50
#> 84: CL00093  cellline_IN tissue_y                    54
#> 85: CL00094  cellline_JO tissue_z                    58
#> 86: CL00095  cellline_KO tissue_z                    62
#> 87: CL00096  cellline_LO tissue_z                    66
#> 88: CL00097  cellline_MO tissue_z                    70
#> 89: CL00098  cellline_NO tissue_z                    74
#> 90: CL00099  cellline_OO tissue_z                    78
#>        clid CellLineName   Tissue ReferenceDivisionTime
#>      <char>       <char>   <char>                 <num>
```
