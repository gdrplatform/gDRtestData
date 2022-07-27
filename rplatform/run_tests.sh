#/bin/sh
repo_path="$1"

echo "Executing $0"
echo "Environment: ${rp_env}"
echo "Working directory: `pwd`"
echo "Working directory contains: `ls | tr '\n' ' '`"

# exit when any command fails
set -e

echo ">>>>>>>> Running linter"
Rscript -e "gDRstyle::lintPkgDirs('/mnt/vol/gDRtestData')"

echo ">>>>> RUNNING UNIT TESTS"
Rscript -e "testthat::test_local(path = '$repo_path', stop_on_failure = TRUE)"

echo ">>>>> RUNNING CHECK"
R CMD build /mnt/vol/gDRtestData &&
    R CMD check gDRtestData_*.tar.gz --no-vignettes --no-examples --no-manual --no-tests

echo ">>>>>>>> RUNNING CHECK DEPENDENCIES"
Rscript -e "gDRstyle::checkDependencies(desc_path='/mnt/vol/gDRtestData/DESCRIPTION', dep_path='/mnt/vol/rplatform/dependencies.yaml')"
