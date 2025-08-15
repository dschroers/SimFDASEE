#' Determine the Number of Principal Components to Reach a Given Variance Threshold
#'
#' This function calculates the minimum number of principal components required
#' to explain at least a specified proportion of the total variance of a covariance matrix.
#'
#' @param C Numeric matrix. Covariance matrix whose eigenvalues are analyzed.
#' @param rho Numeric between 0 and 1. Desired proportion of total variance to explain. Default is 0.9.
#'
#' @return Integer. The minimum number of components needed to reach at least `rho` proportion of total variance.
#'
#' @examples
#' C <- matrix(c(2, 1, 1, 2), 2, 2)
#' d_star(C, rho = 0.9)
#'
#' @export
d_star <- function(C, rho = 0.9) {
  eigvals <- eigen(C, only.values = TRUE)$values
  scores <- cumsum(eigvals) / sum(eigvals)
  return(min(which(scores >= rho)))
}
