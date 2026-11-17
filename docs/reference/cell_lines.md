# Cell lines

Cell lines

## Value

data.table

## Examples

``` r
path <- system.file("annotation_data", "cell_lines.csv", package = "gDRtestData")
data.table::fread(file = path)
#>        clid               CellLineName Tissue ReferenceDivisionTime
#>      <char>                     <char> <char>                 <int>
#>  1: CL00016                cellline_GB breast                    46
#>  2: CL00017                cellline_HB breast                    50
#>  3: CL00011                cellline_BA breast                    26
#>  4: CL00012                cellline_CA breast                    30
#>  5: CL00013                cellline_DA breast                    34
#>  6: CL00014                cellline_EA breast                    38
#>  7: CL00015                cellline_FA breast                    42
#>  8: CL00018                cellline_IB breast                    54
#>  9: CL00019                cellline_JB breast                    58
#> 10: CL00020                cellline_KB breast                    62
#> 11: CL00100                     CAMA-1 Breast                    53
#> 12: CL00101                     EFM-19 Breast                    64
#> 13: CL00102                      MCF-7 Breast                    30
#> 14: CL00103                      T-47D Breast                    60
#> 15: CL00104                       501A   Skin                    27
#> 16: CL00105                    624 mel   Skin                    26
#> 17: CL00106                    928 mel   Skin                    27
#> 18: CL00107                      A-253   Skin                    47
#> 19: CL00108                      A-375   Skin                    25
#> 20: CL00109                      A-431   Skin                    30
#> 21: CL00110                      A2058   Skin                    18
#> 22: CL00111                       A388   Skin                    45
#> 23: CL00112                   COLO 800   Skin                    45
#> 24: CL00113                   COLO 829   Skin                    56
#> 25: CL00114                   COLO 849   Skin                    47
#> 26: CL00115                   COLO 853   Skin                    NA
#> 27: CL00116                   COLO 857   Skin                    31
#> 28: CL00117                   COLO 858   Skin                    NA
#> 29: CL00118                      G-361   Skin                    35
#> 30: CL00119                    HCC1498   Skin                    31
#> 31: CL00120                      HMY-1   Skin                    31
#> 32: CL00121                      HSC-1   Skin                    57
#> 33: CL00122                      HSC-5   Skin                    NA
#> 34: CL00123                     HT-144   Skin                    57
#> 35: CL00124                    Hs 294T   Skin                    36
#> 36: CL00125                    Hs 695T   Skin                    56
#> 37: CL00126                   Hs 852.T   Skin                    63
#> 38: CL00127                      IGR-1   Skin                    18
#> 39: CL00128                     IGR-39   Skin                    45
#> 40: CL00129                 MDA-MB-435   Skin                    35
#> 41: CL00130                     MEL-HO   Skin                    27
#> 42: CL00131                   MEL-JUSO   Skin                    29
#> 43: CL00132                       MMAc   Skin                    82
#> 44: CL00133                  RPMI-7951   Skin                    30
#> 45: CL00134                    RVH-421   Skin                    50
#> 46: CL00135                   SK-MEL-1   Skin                    63
#> 47: CL00136                   SK-MEL-2   Skin                    62
#> 48: CL00137                  SK-MEL-24   Skin                    72
#> 49: CL00138                  SK-MEL-28   Skin                    39
#> 50: CL00139                   SK-MEL-3   Skin                    50
#> 51: CL00140                  SK-MEL-30   Skin                    26
#> 52: CL00141                  SK-MEL-31   Skin                    55
#> 53: CL00142                   SK23-mel   Skin                    20
#> 54: CL00143                   UACC-257   Skin                    47
#> 55: CL00144                    UACC-62   Skin                    39
#> 56: CL00145                  UCSD-242l   Skin                    46
#> 57: CL00146                   WM-266-4   Skin                    29
#> 58: CL00147                    NCI-H23   Lung                    37
#> 59: CL00148                   NCI-H358   Lung                    50
#> 60: CL00149  NCI-H23.AMG-510resistance   Lung                    NA
#> 61: CL00150 NCI-H358.AMG-510resistance   Lung                    NA
#>        clid               CellLineName Tissue ReferenceDivisionTime
#>      <char>                     <char> <char>                 <int>
#>     parental_identifier                 subtype
#>                  <char>                  <char>
#>  1:         cellline_GB                    <NA>
#>  2:         cellline_HB                    <NA>
#>  3:         cellline_BA                    <NA>
#>  4:         cellline_CA                    <NA>
#>  5:         cellline_DA                    <NA>
#>  6:         cellline_EA                    <NA>
#>  7:         cellline_FA                    <NA>
#>  8:         cellline_IB                    <NA>
#>  9:         cellline_JB                    <NA>
#> 10:         cellline_KB                    <NA>
#> 11:              CAMA-1          Adenocarcinoma
#> 12:              EFM-19               Carcinoma
#> 13:               MCF-7          Adenocarcinoma
#> 14:               T-47D   Adenocarcinoma Ductal
#> 15:                501A                Melanoma
#> 16:                501A                Melanoma
#> 17:             928 mel                Melanoma
#> 18:               A-253               Carcinoma
#> 19:               A-375                Melanoma
#> 20:               A-431 Carcinoma Squamous Cell
#> 21:               A2058                Melanoma
#> 22:                A388    Carcinoma Metastatic
#> 23:            COLO 800                Melanoma
#> 24:            COLO 829                Melanoma
#> 25:            COLO 849     Melanoma Metastatic
#> 26:            COLO 853                Melanoma
#> 27:            COLO 857                Melanoma
#> 28:            COLO 858                Melanoma
#> 29:               G-361                Melanoma
#> 30:             HCC1498                Melanoma
#> 31:               HMY-1                Melanoma
#> 32:               HSC-1 Carcinoma Squamous Cell
#> 33:               HSC-5 Carcinoma Squamous Cell
#> 34:              HT-144                Melanoma
#> 35:               A101D                Melanoma
#> 36:             Hs 695T                Melanoma
#> 37:            Hs 852.T                Melanoma
#> 38:               IGR-1                Melanoma
#> 39:              IGR-39                Melanoma
#> 40:          MDA-MB-435                Melanoma
#> 41:              MEL-HO                Melanoma
#> 42:            MEL-JUSO                Melanoma
#> 43:                MMAc                Melanoma
#> 44:           RPMI-7951                Melanoma
#> 45:             RVH-421                Melanoma
#> 46:            SK-MEL-1     Melanoma Metastatic
#> 47:            SK-MEL-2                Melanoma
#> 48:           SK-MEL-24                Melanoma
#> 49:           SK-MEL-28                Melanoma
#> 50:            SK-MEL-3                Melanoma
#> 51:           SK-MEL-30                Melanoma
#> 52:           SK-MEL-31                Melanoma
#> 53:           SK-23 MEL                Melanoma
#> 54:            UACC-257                Melanoma
#> 55:             UACC-62     Melanoma Metastatic
#> 56:           UCSD-242l                Melanoma
#> 57:            WM-266-4                Melanoma
#> 58:             NCI-H23                   NSCLC
#> 59:            NCI-H358                   NSCLC
#> 60:             NCI-H23                   NSCLC
#> 61:            NCI-H358                   NSCLC
#>     parental_identifier                 subtype
#>                  <char>                  <char>
```
