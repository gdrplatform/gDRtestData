# get_test_dataset_paths

Returns named vector of absolute paths to test datasets.

## Usage

``` r
get_test_dataset_paths(datasets_dir = NULL, pattern = "finalMAE_")
```

## Arguments

- datasets_dir:

  path to directory with datasets (default `NULL`). If `NULL`, then
  `inst/testdata` directory from `gDRtestData` will be used.

- pattern:

  used to: (1) filter to qs2 files from the dataset_dir path and (2)
  prettify the labels of the files

## Value

named vector of absolute paths

## Author

Kamil Foltyński <kamil.foltynski@contractors.roche.com>

## Examples

``` r
get_test_dataset_paths()
#>                                                                                 combo_2dose_nonoise 
#>    "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_2dose_nonoise.qs2" 
#>                                                                                combo_2dose_nonoise2 
#>   "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_2dose_nonoise2.qs2" 
#>                                                                                combo_2dose_nonoise3 
#>   "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_2dose_nonoise3.qs2" 
#>                                                                                    combo_codilution 
#>       "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_codilution.qs2" 
#>                                                                              combo_codilution_small 
#> "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_codilution_small.qs2" 
#>                                                                                        combo_matrix 
#>           "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_matrix.qs2" 
#>                                                                                  combo_matrix_small 
#>     "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_matrix_small.qs2" 
#>                                                                                        combo_triple 
#>           "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_triple.qs2" 
#>                                                                                              medium 
#>                 "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_medium.qs2" 
#>                                                                                               small 
#>                  "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_small.qs2" 
#>                                                                                      small_no_noise 
#>         "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_small_no_noise.qs2" 
#>                                                                                             wLigand 
#>                "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_wLigand.qs2" 
path <- system.file("testdata", package = "gDRtestData", mustWork = TRUE)
get_test_dataset_paths(path)
#>                                                                                 combo_2dose_nonoise 
#>    "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_2dose_nonoise.qs2" 
#>                                                                                combo_2dose_nonoise2 
#>   "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_2dose_nonoise2.qs2" 
#>                                                                                combo_2dose_nonoise3 
#>   "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_2dose_nonoise3.qs2" 
#>                                                                                    combo_codilution 
#>       "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_codilution.qs2" 
#>                                                                              combo_codilution_small 
#> "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_codilution_small.qs2" 
#>                                                                                        combo_matrix 
#>           "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_matrix.qs2" 
#>                                                                                  combo_matrix_small 
#>     "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_matrix_small.qs2" 
#>                                                                                        combo_triple 
#>           "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_combo_triple.qs2" 
#>                                                                                              medium 
#>                 "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_medium.qs2" 
#>                                                                                               small 
#>                  "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_small.qs2" 
#>                                                                                      small_no_noise 
#>         "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_small_no_noise.qs2" 
#>                                                                                             wLigand 
#>                "/tmp/RtmpvM8NKN/temp_libpath21564ef55047/gDRtestData/testdata/finalMAE_wLigand.qs2" 
```
