# gDRtestData

## Overview

The `gDRtestData` package is intended to store and generate example data
that can be used through the `gDR` suite.

## Use Cases

### Synthetic data generation

Since gDR is a computational suite for drug response data from any
experiment, a synthetic dataset is also needed for testing and
exploration.

The basis of this package are two functions to generate the synthetic
sets of cell lines and drugs.

``` r
cell_lines <- create_synthetic_cell_lines()
drugs <- create_synthetic_drugs()
```

These base objects can be extended with additional information.

1.  Replicates

``` r
cl_rep <- add_data_replicates(cell_lines)
head(cl_rep)
```

2.  Concentration

``` r
cl_conc <- add_concentration(cell_lines)
head(cl_conc)
```

Or the user can do both with one function:

``` r
df_layout <- prepareData(cell_lines, drugs)
head(df_layout)
```

Additionally, the user may fill in the full response data with the day 0
part.

``` r
df_layout_small <- prepareData(cell_lines[seq_len(2), ], drugs[seq_len(4), ])
df_layout_small$Duration <- 72
df_layout_small$ReadoutValue <- 0
df_layout_small_with_Day0 <- add_day0_data(df_layout_small)
head(df_layout_small_with_Day0)
```

In a further step, the user may generate a set of synthetic results:

1.  Hill coefficient

``` r
hill <- generate_hill_coef(cell_lines, drugs)
```

2.  EC50 metric

``` r
ec50_met <- generate_ec50(cell_lines, drugs)
```

3.  E inf metric

``` r
einf_met <- generate_e_inf(cell_lines, drugs)
```

Or the user can create full response data with one function (for
single-agent):

``` r
response_data <- prepareMergedData(cell_lines, drugs)
head(response_data)
```

**SUMMARY**

| Step | Function                      | Output (`data.table`)                                        |
|:----:|:------------------------------|:-------------------------------------------------------------|
|  0   | create_synthetic_cell_lines() | synthetic cell lines                                         |
|  0   | create_synthetic_drugs()      | synthetic drugs                                              |
|  1   | prepareData()                 | cell lines and drug merged with replicates and concentration |
|  2   | prepareMergedData()           | full response data for single-agent                          |
|  2   | prepareComboMergedData()      | full response data for combo                                 |
|  2   | prepareCodilutionData ()      | full response data for co-dilution                           |

### Synthetic object of gDR data model

The gDR data model is built on the MultiAssayExperiments (MAE)
structure. A detailed description of the gDR data model can be found in
`gDRcore` package vignette.

In `inst/testdata` the user may find a set of `qs2` files that are
examples of gDR data model for different data types. In the file
`synthetic_list.yml` one can find a list of these datasets. Currently
available are:

- combo_2dose_nonoise,
- combo_2dose_nonoise2,
- combo_2dose_nonoise3,
- combo_codilution_small,
- combo_codilution,
- combo_matrix_small,
- combo_matrix,
- combo_triple,
- medium,
- small_no_noise,
- small,
- wLigand,
- prism

The script `generate_example_data.R` which shows how to generate and
process above-mentioned datasets is in `inst/scripts` dir. All key
functions can be found in package `gDRcore` in script
`generate_wrappers.R`.

*Note*: PRISM data is not created synthetically by the script.

### Annotation data

In `inst/annotation_data` the user can find `CSV` files used in
`gDRcore` for testing annotation functions.

### Other

Other files which were not mentioned above are used for testing gDR
suite functionality.

## SessionInfo

``` r
sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] gDRtestData_1.11.6 BiocStyle_2.40.0  
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.6           knitr_1.51          rlang_1.2.0        
#>  [4] xfun_0.59           otel_0.2.0          textshaping_1.0.5  
#>  [7] data.table_1.18.4   jsonlite_2.0.0      backports_1.5.1    
#> [10] htmltools_0.5.9     ragg_1.5.2          sass_0.4.10        
#> [13] rmarkdown_2.31      evaluate_1.0.5      jquerylib_0.1.4    
#> [16] fastmap_1.2.0       yaml_2.3.12         lifecycle_1.0.5    
#> [19] bookdown_0.47       BiocManager_1.30.27 compiler_4.6.0     
#> [22] fs_2.1.0            systemfonts_1.3.2   digest_0.6.39      
#> [25] R6_2.6.1            checkmate_2.3.4     bslib_0.11.0       
#> [28] tools_4.6.0         pkgdown_2.2.0       cachem_1.1.0       
#> [31] desc_1.4.3
```
