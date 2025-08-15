#' Functional Data Truncation Based on Covariance Eigenstructure
#'
#' This function computes a truncation function for functional data and identifies
#' time points where the function exceeds a threshold based on the covariance structure.
#' It is useful for detecting extreme or influential observations in functional time series.
#'
#' @param d Integer. Number of leading eigenvalues/components to consider for inversion.
#' @param C Numeric matrix. Covariance matrix of the functional data.
#' @param data Numeric matrix. Functional data with rows corresponding to time points and columns to spatial/grid points.
#' @param Delta Numeric. Increment size used in the threshold calculation.
#' @param sd Numeric. Number of standard deviations used to define the truncation threshold.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{gn}{Numeric vector of truncation function values for each time point.}
#'   \item{locations}{Integer indices of time points where `gn` exceeds the threshold.}
#' }
#'
#' @examples
#' x_grid <- seq(0, 1, length.out = 100)
#' q <- outer(x_grid, x_grid, function(x, y) exp(-abs(x - y)))
#' data <- matrix(rnorm(100*50), nrow = 50, ncol = 100)
#' functional_data_truncation(d = 5, C = q, data = data, Delta = 0.01, sd = 3)
#'
#' @export
functional_data_truncation <- function(d, C, data, Delta, sd) {
  E <- eigen(C)
  n <- nrow(data)

  # Calculate the values of the truncation function
  gn <- numeric(n)
  for (i in 1:n) {
    x <- data[i, ]

    # First part of g_n (sum of first 'd' components)
    g1 <- sum((x %*% E$vectors[, 1:d])^2 / E$values[1:d])

    # Second part of g_n (remaining components, if any)
    idx <- seq(d + 1, length(E$values))
    g2_num <- sum((x %*% E$vectors[, idx])^2)
    g2_den <- sum(E$values[idx])
    g2 <- g2_num / g2_den

    # Take the square root of their sum
    gn[i] <- sqrt(g1 + g2)
  }

  # Identify truncation locations
  truncation.locations <- which(gn > (sd * sqrt(d + 1) * Delta^(0.49)))

  return(list("gn" = gn, "locations" = truncation.locations))
}
