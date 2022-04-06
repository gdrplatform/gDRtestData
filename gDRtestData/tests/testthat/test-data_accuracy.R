library(SummarizedExperiment)
library(BumpyMatrix)
library(gDRutils)
library(gDRwrapper)
library(gDRcore)
library(reshape2)

cell_lines <- create_synthetic_cell_lines()
drugs <- create_synthetic_drugs()
e_inf <- generate_e_inf(drugs, cell_lines)
ec50 <- generate_ec50(drugs, cell_lines)
hill_coef <- generate_hill_coef(drugs, cell_lines)

evaluateData <- function(data, e_inf, ec50, hill_coef, vals, cols = NULL, FUN = all) {
  dt_test <- test_accuracy(data, e_inf, ec50, hill_coef)
  
  dt_test <- ifelse(is.null(cols), dt_test, dt_test[cols, ])
  
  apply(
    abs(dt_test) < vals, 1, 
    function(x) testthat::expect_true(FUN(x))
  )
}

evaluateMeta <- function(data, expectedLength, type = "fixed") {
  combos <- metadata(data)$drug_combinations
  testthat::expect_equal(length(combos), expectedLength)
  testthat::expect_true(
    all(sapply(
      combos, 
      function(x) length(x$condition$Concentration_2) == (length(x$rows) - 1)
    ))
  )
  testthat::expect_true(all(sapply(combos, "[[", "type") == type))
}

evaluateComboTable <- function(data, x1, x2, x3 = NULL) {
  rows <- rowData(data[[1]])
  testthat::expect_true(
    all(table(rows[, c("Gnumber", "Gnumber_2")])[, drugs$Gnumber[c(21, 26)]] == x1)
  )
  testthat::expect_true(
    all(table(rows[rows$DrugName_2 != "vehicle","Concentration_2"]) == x2)
  )
  if (!is.null(x3)) {
    testthat::expect_true(
      all(table(rows[, paste0("Concentration", c("_2", "_3"))]) == x3)
    )
  }
}

evaluateComboDt <- function(data, x1, x2) {
  dt <- convert_mae_assay_to_dt(data, "Averaged")
  table <- table(dt[dt$DrugName_2 != "vehicle", c("Concentration", "Concentration_2")])
  testthat::expect_true(all(dim(table) == x1))
  testthat::expect_true(all(table == x2))
}

getDelta <- function(mae1, mae2, cols) {
  DT1 <- convert_mae_assay_to_dt(mae1, "Metrics")
  DT2 <- convert_mae_assay_to_dt(mae2, "Metrics")
  DT1$rId <- gsub("_\\d\\d?$", "", DT1$rId)
  DT2$rId <- gsub("_\\d\\d?$", "", DT2$rId)
  
  merge(
    DT1[ , c("rId", "cId", "normalization_type", "x_0", cols)],
    DT2[ , c("rId", "cId", "normalization_type", "x_0", cols)], 
    by = c("rId", "cId", "normalization_type")
  )
}

testthat::context("Testing generated data accuracy")

testthat::test_that(
  desc = "No noise raw data",
  code = {
    # test accurarcy of the processing and fitting (no noise => low tolerance)
    data <- generateNoNoiseRawData(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
    evaluateData(data, e_inf, ec50, hill_coef, c(1e-3, 2.2e-3, 0.04, 0.015, 1e-4))
  }
)

testthat::test_that(
  desc = "Noise raw data",
  code = {
    # test accurarcy of the processing and fitting (noise => medium tolerance)
    data <- generateNoiseRawData(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
    evaluateData(data, e_inf, ec50, hill_coef, c(0.5, 0.1, 1.5, 1.2, 0.05))
  }
)

testthat::test_that(
  desc = "Ligand data",
  code = {
    # test accurarcy of the processing and fitting for Ligand = 0.1 (noise => medium tolerance)
    data <- generateLigandData(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
    rows <- SummarizedExperiment::rowData(data)

    evaluateData(
      data[rows$Ligand > 0, ], 
      e_inf, 
      ec50, 
      hill_coef, 
      c(1e-3, 3e-3, 0.031, 0.015, 1e-4)
    )
    # test fit quality for Ligand = 0 and that delta(e_inf) < 0
    evaluateData(
      data[rows$Ligand == 0, ],
      e_inf,
      ec50,
      hill_coef,
      vals = c(-0.15, 1e-4),
      cols = c("delta_einf", "1_r2")
    )
  }
)

testthat::test_that(
  desc = "Medium data",
  code = {
    # test accurarcy of the processing and fitting (noise => medium tolerance)
    data <- generateMediumData(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
    evaluateData(data, e_inf, ec50, hill_coef, c(0.5, 0.2, 2.5, 1.2, 0.3))
  }
)

testthat::test_that(
  desc = "Many lines data",
  code = {
    # test accurarcy of the processing and fitting (noise => medium tolerance)
    data <- generateManyLinesData(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
    evaluateData(data, e_inf, ec50, hill_coef, c(0.5, 0.2, 2.5, 1.2, 0.3))
  }
)

testthat::test_that(
  desc = "Many drugs data",
  code = {
    # test accurarcy of the processing and fitting (noise => medium tolerance)
    data <- generateManyDrugsData(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
    evaluateData(data, e_inf, ec50, hill_coef, c(0.5, 0.2, 2.5, 1.2, 0.3))
  }
)

testthat::test_that(
  desc = "Combo no noise data",
  code = {
    ## 1st case
    data <- generateComboNoNoiseData(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
    evaluateData(data[[1]], e_inf, ec50, hill_coef, c(1e-3, 2e-3, 0.02, 0.015, 1e-4))
    evaluateMeta(data[[1]], 3)
    
    ## 2nd case
    data2 <- generateComboNoNoiseData2(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
    evaluateData(
      data2[[1]][rowData(data2[[1]])$Concentration_2 == 0, ], 
      e_inf, 
      ec50, 
      hill_coef, 
      c(1e-3, 2e-3, 0.02, 0.015, 1e-4)
    )

    # compare to other way of processing the data
    delta <- getDelta(data, data2, c("x_max"))
    testthat::expect_true(
      all(abs(quantile((delta$x_0.x - delta$x_0.y)[!grepl("vehicle", delta$rId)])) < .0005)
    )
    
    ## 3rd case
    data3 <- generateComboNoNoiseData3(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
    evaluateData(
      data3[[1]][rowData(data3[[1]])$Concentration_2 == 0, ], 
      e_inf, 
      ec50, 
      hill_coef, 
      c(1e-3, 2e-3, 0.02, 0.015, 1e-4)
    )
    
    # compare to the complete data
    delta2 <- getDelta(data, data3, c("x_inf", "r2"))
    # checking the x_0 value was properly fitted for most cases (there are a few failures but it is ok)
    testthat::expect_lt(sum(abs(delta2$x_0.x - delta2$x_0.y) > .004), 4)
    testthat::expect_lt(sum(abs(delta2$x_0.x - delta2$x_0.y) > .025), 3)
    # checking the x_inf value was properly fitted
    testthat::expect_lt(sum(abs(delta2$x_inf.x - delta2$x_inf.y) > .00003), 5)
    testthat::expect_lt(sum(abs(delta2$x_inf.x - delta2$x_inf.y) > .008), 3)
  }
)

testthat::test_that(
  desc = "Combo many drugs data",
  code = {
    # test accuracy of the processing and fitting for the single agent
    data <- generateComboManyDrugs(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
    evaluateData(
      data[[1]][rowData(mae[[1]])$Concentration_2 == 0, ], 
      e_inf, 
      ec50, 
      hill_coef, 
      vals = c(0.5, 0.2, 2.5, 1.2, 0.3)
    )
    # test the effect of the combination treatment
    evalFun <- function(x) {
      sum(x) == 2
    }
    evaluateData(
      data[[1]][rowData(mae[[1]])$Concentration_2 > 0, ], 
      e_inf, 
      ec50, 
      hill_coef, 
      vals = c(-.1, .01),
      cols = c("delta_einf", "1_r2"),
      FUN = evalFun
    )
    evaluateMeta(data[[1]], 149)
  }
)

testthat::test_that(
  desc = "Combo matrix data",
  code = {
    ## Small matrix
    # test accuracy of the processing and fitting for the single agent
    data <- generateComboMatrixSmall(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
    evaluateData(
      data[[1]][rowData(data[[1]])$Concentration_2 == 0, ],
      e_inf, 
      ec50, 
      hill_coef, 
      c(1e-3, 6e-3, 0.12, 0.015, 1e-4)
    )
    evaluateComboTable(data, 8, 6)
    evaluateComboDt(data, 8, 36)
    evaluateMeta(data[[1]], 6, "matrix")
    # add and test calculation for combo matrix
    # TODO when the functions are cleaned up
    
    ## Mid-size matrix
    data2 <- generateComboMatrix(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
    evaluateComboTable(data2, 9, 18)
    evaluateComboDt(data2, 9, 144)
    evaluateMeta(data2[[1]], 18, "matrix")
    # add and test calculation for combo matrix
    # TODO when the functions are cleaned up
  }
)

testthat::test_that(
  desc = "Triple combo data",
  code = {
    # test accuracy of the processing and fitting for the single agent
    data <- generateTripleComboMatrix(cell_lines, drugs, e_inf, ec50, hill_coef, FALSE)
    
    dt_test <- data[[1]][rowData(data[[1]])$Concentration_2 == 0 &
      rowData(data[[1]])$Concentration_3 == 0, ]
    evaluateData(dt_test, e_inf, ec50, hill_coef, c(1e-3, 6e-3, 0.12, 0.015, 1e-4))
    evaluateComboTable(data, 24, 18, c(3, array(6, 8)))

    dt <- convert_mae_assay_to_dt(data, "Averaged")
    table <- table(dt[dt$DrugName_2 != "vehicle", paste0("Concentration",c("", "_2", "_3"))])
    testthat::expect_true(all(dim(table) == c(8, 8, 3)))
    
    evaluateMeta(data[[1]], 18, "matrix")
    # add and test calculation for combo matrix
    # TODO when the functions are cleaned up
  }
)
