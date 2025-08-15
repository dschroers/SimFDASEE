#' Quantile-Based Truncation of Functional Data
#'
#' This function identifies observations (rows) in a functional data matrix
#' whose Euclidean norm exceeds a specified quantile threshold. These observations
#' are considered extreme and can be removed for robust statistical analysis.
#'
#' @param data Numeric matrix. Functional data with rows corresponding to time points and columns to spatial/grid points.
#' @param q Numeric between 0 and 1. Quantile threshold; rows with norms greater than or equal to this quantile are truncated.
#'
#' @return Integer vector of row indices of `data` exceeding the specified quantile threshold.
#'
#' @details
#' The Euclidean norm of each row is computed to measure the magnitude of the functional observation.
#' Rows with norms greater than or equal to the `q`-th quantile are returned as truncation locations.
#' This is used as a preprocessing step for robust covariance estimation or functional data truncation.
#'
#' @examples
#' # Simulate functional data
#' set.seed(123)
#' data <- matrix(rnorm(50*100), nrow = 50, ncol = 100)
#' quantile_truncation(data, q = 0.95)
#'
#' @export
quantile_truncation <- function(data, q) {

  # Compute Euclidean norm for each row
  norms.data <- apply(data, 1, function(row) sqrt(sum(row^2)))

  # Determine quantile threshold
  quant <- quantile(norms.data, probs = q)

  # Identify indices of points exceeding the threshold
  truncation.locations <- which(norms.data >= quant)
  return(truncation.locations)
}
