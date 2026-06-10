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
#> 11: CL00100       CAMA-1 Breast                    53              CAMA-1
#> 12: CL00101       EFM-19 Breast                    64              EFM-19
#> 13: CL00102        MCF-7 Breast                    30               MCF-7
#> 14: CL00103        T-47D Breast                    60               T-47D
#> 15: CL00104         501A   Skin                    27                501A
#> 16: CL00105      624 mel   Skin                    26                501A
#> 17: CL00106      928 mel   Skin                    27             928 mel
#> 18: CL00107        A-253   Skin                    47               A-253
#> 19: CL00108        A-375   Skin                    25               A-375
#> 20: CL00109        A-431   Skin                    30               A-431
#> 21: CL00110        A2058   Skin                    18               A2058
#> 22: CL00111         A388   Skin                    45                A388
#> 23: CL00112     COLO 800   Skin                    45            COLO 800
#> 24: CL00113     COLO 829   Skin                    56            COLO 829
#> 25: CL00114     COLO 849   Skin                    47            COLO 849
#> 26: CL00115     COLO 853   Skin                    NA            COLO 853
#> 27: CL00116     COLO 857   Skin                    31            COLO 857
#> 28: CL00117     COLO 858   Skin                    NA            COLO 858
#> 29: CL00118        G-361   Skin                    35               G-361
#> 30: CL00119      HCC1498   Skin                    31             HCC1498
#> 31: CL00120        HMY-1   Skin                    31               HMY-1
#> 32: CL00121        HSC-1   Skin                    57               HSC-1
#> 33: CL00122        HSC-5   Skin                    NA               HSC-5
#> 34: CL00123       HT-144   Skin                    57              HT-144
#> 35: CL00124      Hs 294T   Skin                    36               A101D
#> 36: CL00125      Hs 695T   Skin                    56             Hs 695T
#> 37: CL00126     Hs 852.T   Skin                    63            Hs 852.T
#> 38: CL00127        IGR-1   Skin                    18               IGR-1
#> 39: CL00128       IGR-39   Skin                    45              IGR-39
#> 40: CL00129   MDA-MB-435   Skin                    35          MDA-MB-435
#> 41: CL00130       MEL-HO   Skin                    27              MEL-HO
#> 42: CL00131     MEL-JUSO   Skin                    29            MEL-JUSO
#> 43: CL00132         MMAc   Skin                    82                MMAc
#> 44: CL00133    RPMI-7951   Skin                    30           RPMI-7951
#> 45: CL00134      RVH-421   Skin                    50             RVH-421
#> 46: CL00135     SK-MEL-1   Skin                    63            SK-MEL-1
#> 47: CL00136     SK-MEL-2   Skin                    62            SK-MEL-2
#> 48: CL00137    SK-MEL-24   Skin                    72           SK-MEL-24
#> 49: CL00138    SK-MEL-28   Skin                    39           SK-MEL-28
#> 50: CL00139     SK-MEL-3   Skin                    50            SK-MEL-3
#> 51: CL00140    SK-MEL-30   Skin                    26           SK-MEL-30
#> 52: CL00141    SK-MEL-31   Skin                    55           SK-MEL-31
#> 53: CL00142     SK23-mel   Skin                    20           SK-23 MEL
#> 54: CL00143     UACC-257   Skin                    47            UACC-257
#> 55: CL00144      UACC-62   Skin                    39             UACC-62
#> 56: CL00145    UCSD-242l   Skin                    46           UCSD-242l
#> 57: CL00146     WM-266-4   Skin                    29            WM-266-4
#>        clid CellLineName Tissue ReferenceDivisionTime parental_identifier
#>      <char>       <char> <char>                 <int>              <char>
#>                     subtype
#>                      <char>
#>  1:                    <NA>
#>  2:                    <NA>
#>  3:                    <NA>
#>  4:                    <NA>
#>  5:                    <NA>
#>  6:                    <NA>
#>  7:                    <NA>
#>  8:                    <NA>
#>  9:                    <NA>
#> 10:                    <NA>
#> 11:          Adenocarcinoma
#> 12:               Carcinoma
#> 13:          Adenocarcinoma
#> 14:   Adenocarcinoma Ductal
#> 15:                Melanoma
#> 16:                Melanoma
#> 17:                Melanoma
#> 18:               Carcinoma
#> 19:                Melanoma
#> 20: Carcinoma Squamous Cell
#> 21:                Melanoma
#> 22:    Carcinoma Metastatic
#> 23:                Melanoma
#> 24:                Melanoma
#> 25:     Melanoma Metastatic
#> 26:                Melanoma
#> 27:                Melanoma
#> 28:                Melanoma
#> 29:                Melanoma
#> 30:                Melanoma
#> 31:                Melanoma
#> 32: Carcinoma Squamous Cell
#> 33: Carcinoma Squamous Cell
#> 34:                Melanoma
#> 35:                Melanoma
#> 36:                Melanoma
#> 37:                Melanoma
#> 38:                Melanoma
#> 39:                Melanoma
#> 40:                Melanoma
#> 41:                Melanoma
#> 42:                Melanoma
#> 43:                Melanoma
#> 44:                Melanoma
#> 45:                Melanoma
#> 46:     Melanoma Metastatic
#> 47:                Melanoma
#> 48:                Melanoma
#> 49:                Melanoma
#> 50:                Melanoma
#> 51:                Melanoma
#> 52:                Melanoma
#> 53:                Melanoma
#> 54:                Melanoma
#> 55:     Melanoma Metastatic
#> 56:                Melanoma
#> 57:                Melanoma
#>                     subtype
#>                      <char>
```
