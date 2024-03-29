% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cavi_loading.R
\name{compute_loading_row_conditional_mean_weight_matrix}
\alias{compute_loading_row_conditional_mean_weight_matrix}
\title{Loading Row Conditional Mean Weight}
\usage{
compute_loading_row_conditional_mean_weight_matrix(
  loading_prior_correlation,
  perform_checks = TRUE
)
}
\arguments{
\item{loading_prior_correlation}{A D'xD' correlation matrix. The correlation
within loadings under the Gaussian prior.}

\item{perform_checks}{Logical. Check if function inputs are specified
correctly.}
}
\value{
A D'x(D'-1) matrix of positive real numbers. The weight of the
  remaining loading matrix in the conditional mean of each row in the loading
  matrix.
}
\description{
Compute the weight of the remaining loading matrix in the conditional mean
for each row of the loading matrix.
}
