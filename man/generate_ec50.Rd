% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random_data.R
\name{generate_ec50}
\alias{generate_ec50}
\title{Calculate EC50 metric}
\usage{
generate_ec50(drugs, cell_lines)
}
\arguments{
\item{drugs}{data.table with drugs}

\item{cell_lines}{data.table with cell lines}
}
\value{
matrix with random EC50
}
\description{
Calculate EC50 metric
}
\examples{
generate_ec50(create_synthetic_drugs(), create_synthetic_cell_lines())

}
\keyword{generate_results}
