#' Preliminary Covariation Estimator for Functional Data
#'
#' This function computes a preliminary estimator of the covariance (or covariation)
#' matrix for functional data after removing extreme observations based on a
#' truncation threshold. It also scales the covariance matrix using a robust
#' estimate of the leading eigenvalue contribution.
#'
#' @param data Numeric matrix. Functional data with rows corresponding to time points and columns to spatial/grid points.
#' @param tq Numeric. Truncation quantile used to identify extreme observations via `quantile_truncation`.
#'
#' @return A numeric matrix representing the preliminary estimated covariance matrix.
#'
#' @details
#' The function first identifies "rough" truncation locations based on the specified quantile `tq`
#' using the `quantile_truncation` function. These rows are removed from the data to reduce the
#' effect of extreme observations. The covariance of the truncated data is computed and scaled
#' by `rho.star`, which is derived from the interquartile range of projections onto the leading
#' eigenvector, to produce a robust preliminary covariation estimator.
#'
#' @examples
#' # Simulate functional data
#' set.seed(123)
#' data <- matrix(rnorm(50*100), nrow = 50, ncol = 100)
#' preliminary_covariation_estimator(data, tq = 0.95)
#'
#' @export
preliminary_covariation_estimator <- function(data, tq) {
  n <- nrow(data)

  # Identify rough truncation locations
  rough.locs <- quantile_truncation(data, tq)
  data_rough_trunc <- data[-rough.locs, , drop = FALSE]

  # Compute preliminary covariance
  C.Prel <- t(data_rough_trunc) %*% data_rough_trunc
  EG <- eigen(C.Prel)

  # Compute robust scaling factor rho.star
  values <- sapply(1:(n - 1), function(i) t(EG$vectors[, 1]) %*% data[i, ])
  q.75 <- quantile(values, 0.75)
  q.25 <- quantile(values, 0.25)

  rho.star <- ((q.75 - q.25)^2) / ((4 * (qnorm(0.75)^2) * EG$values[1]) / n)

  # Scale the covariance
  C.Prel <- rho.star * C.Prel

  return(C.Prel = C.Prel)
}
