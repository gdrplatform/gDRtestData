#' Create data.frame with synthetic cell lines
#'
#' @return data.frame with synthetic cell lines
#'
#' @examples
#' create_synthetic_cell_lines()
#' 
#' @export
create_synthetic_cell_lines <- function() {
  cell_lines <- data.frame(
    clid = paste0("CL000", 11:25),
    CellLineName = paste0("cellline_", LETTERS[seq_len(15)]),
    Tissue = sort(paste0("tissue_", array(letters[c(24, 24:26, 26)], 15))),
    ReferenceDivisionTime = seq(22, 80, 4)
  )
  cell_lines <- Reduce(rbind, list(cell_lines)[rep(1, times = 6)])
  cell_lines$clid <- paste0("CL000", 9 + (seq_len(nrow(cell_lines))))
  cell_lines$CellLineName <- paste0(cell_lines$CellLineName, sort(array(LETTERS[seq_len(15)], 90)))
  cell_lines$Tissue[16:40] <- "tissue_w"
  cell_lines$Tissue[41:50] <- "tissue_v"
  cell_lines
}

#' Create data.frame with synthetic drugs
#'
#' @return data.frame with synthetic drugs
#' @examples
#' create_synthetic_drugs()
#' 
#' @export
create_synthetic_drugs <- function() {
  drugs <- data.frame(
    Gnumber = paste0("G00", 11:50),
    DrugName = paste0("drug_", 11:50),
    drug_moa = sort(paste0("moa_", array(LETTERS[c(1, seq_len(6), 6)], 40)))
  )
  drugs <- Reduce(rbind, list(drugs)[rep(1, times = 6)])
  drugs$Gnumber <- sprintf("G00%03i", seq_len(nrow(drugs)))
  drugs$DrugName <- sprintf("drug_%03i", seq_len(nrow(drugs)))
  drugs$drug_moa[-seq_len(80)] <- sort(paste0("moa_", array(LETTERS[seq_len(24)], 160)))
  drugs
}


#' Generate hill coefficient
#'
#' @param drugs data.frame with drugs
#' @param cell_lines data.frame with cell lines
#'
#' @return matrix with random hill coefficient
#' @examples
#' generate_hill_coef(create_synthetic_drugs(), create_synthetic_cell_lines())
#' 
#' @export
generate_hill_coef <- function(drugs, cell_lines) {
  hill_coef <- matrix(1.8 + stats::runif(nrow(drugs) * nrow(cell_lines)), nrow(drugs), nrow(cell_lines))
  colnames(hill_coef) <- cell_lines$clid
  rownames(hill_coef) <- drugs$Gnumber
  hill_coef
}

#' Calculate EC50 metric
#'
#' @param drugs data.frame with drugs
#' @param cell_lines data.frame with cell lines
#'
#' @return matrix with random EC50
#' @examples
#' generate_ec50(create_synthetic_drugs(), create_synthetic_cell_lines())
#' 
#' @export
generate_ec50 <- function(drugs, cell_lines) {
  nDrugs <- nrow(drugs)
  nCells <- nrow(cell_lines)
  ec50 <- matrix(stats::runif(nDrugs * nCells) - 0.5, nDrugs, nCells) +
    matrix(sort(rep(seq(-1.2, 0, 0.3), 8)), nDrugs, nCells) +
    t(matrix(seq(-0.4, 0, 0.1), nCells, nDrugs))

  tissue_x <- cell_lines$Tissue == "tissue_x"
  tissue_y <- cell_lines$Tissue == "tissue_y"
  tissue_z <- cell_lines$Tissue == "tissue_z"
  moa_AB <- drugs$drug_moa %in% c("moa_A", "moa_B")
  moa_CD <- drugs$drug_moa %in% c("moa_C", "moa_D")
  moa_E <- drugs$drug_moa %in% "moa_E"
  ec50[, tissue_x] <- ec50[, tissue_x] - 0.5 - 0.3 * stats::runif(sum(tissue_x))
  ec50[moa_AB, tissue_y] <- ec50[moa_AB, tissue_y] - 0.4 * stats::runif(sum(tissue_y)) - 0.4
  ec50[moa_CD, tissue_z] <- ec50[moa_CD, tissue_z] + 0.5 * stats::runif(sum(tissue_z)) - 1
  ec50[moa_E, ] <- ec50[moa_E, ] + 0.6 + 0.5 * stats::runif(sum(moa_E))
  ec50 <- 10 ^ ec50

  colnames(ec50) <- cell_lines$clid
  rownames(ec50) <- drugs$Gnumber
  ec50
}


#' Calculate E inf metric
#'
#' @param drugs data.frame with drugs
#' @param cell_lines data.frame with cell lines
#'
#' @return matrix with random E inf
#' @examples
#' generate_e_inf(create_synthetic_drugs(), create_synthetic_cell_lines())
#' 
#' @export
generate_e_inf <- function(drugs, cell_lines) {
  nDrugs <- nrow(drugs)
  nCells <- nrow(cell_lines)
  
  e_inf <- matrix(0.5 * stats::runif(nDrugs * nCells), nDrugs, nCells) +
    t(matrix(seq(0, 0.2, 0.05), nCells, nDrugs))

  tissue_x <- cell_lines$Tissue == "tissue_x"
  tissue_y <- cell_lines$Tissue == "tissue_y"
  tissue_z <- cell_lines$Tissue == "tissue_z"
  moa_AC <- drugs$drug_moa %in% c("moa_A", "moa_C")
  moa_CE <- drugs$drug_moa %in% c("moa_C", "moa_E")
  moa_BE <- drugs$drug_moa %in% c("moa_B", "moa_E")
  moa_F <- drugs$drug_moa %in% "moa_F"
  e_inf[, tissue_x] <- e_inf[, tissue_x] + 0.3 * stats::runif(sum(tissue_x)) + 0.1
  e_inf[moa_AC, tissue_y] <- e_inf[moa_AC, tissue_y] - 0.2 * stats::runif(sum(tissue_y)) - 0.2
  e_inf[moa_CE, tissue_z] <- e_inf[moa_CE, tissue_z] - 0.5 * stats::runif(sum(tissue_z))
  e_inf[moa_BE, tissue_x] <- e_inf[moa_BE, tissue_x] - 0.5 * stats::runif(sum(tissue_x))
  e_inf[moa_F, ] <- e_inf[moa_F, ] + 0.3 + 0.2 * stats::runif(sum(moa_F))
  e_inf <- matrix(pmin(0.89, pmax(0.01, e_inf)) + stats::runif(nDrugs * nCells) * 0.1, nDrugs, nCells)
  
  colnames(e_inf) <- cell_lines$clid
  rownames(e_inf) <- drugs$Gnumber
  
  e_inf
}
