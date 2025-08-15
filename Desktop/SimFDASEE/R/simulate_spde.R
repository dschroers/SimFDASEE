#' Simulate BMTSM Functional Data with Optional Jumps
#'
#' Generates a matrix of functional data over space and time, combining a continuous
#' part generated from a stationary covariance kernel and an optional discontinuous
#' jump component.
#'
#' @param nx Integer. Number of spatial grid points. Default is 101.
#' @param nt Integer. Number of time points. Default is 101.
#' @param sigma Numeric. Scaling factor for the covariance kernel of the continuous part. Default is 1.
#' @param lambda Numeric. Rate parameter for the exponential distribution controlling jump times. Default is 2.
#' @param kappa Numeric. Variance parameter for jump magnitudes. Default is 0.5.
#' @param kernel Function. Stationary covariance function for the continuous part. Should accept a single argument \(h = x - y\). Default is exponential \code{function(h) exp(-abs(h))}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{cont_samples}{Matrix of continuous functional samples (nt × nx).}
#'   \item{discont_samples}{Matrix including the jump component (nt × nx).}
#'   \item{locs}{Vector of time indices where jumps occur.}
#' }
#'
#' @details
#' The continuous part of the data is generated via a stationary covariance kernel.
#' The jump component is generated at random times from a Poisson-like process,
#' with jump magnitudes drawn from a normal distribution with variance \code{kappa}.
#' The covariance function can be customized to simulate different correlation structures.
#'
#' @examples
#' # Default simulation with exponential kernel
#' sim <- simulate_spde()
#'
#' # Simulation with Gaussian kernel
#' sim2 <- simulate_spde(kernel = function(h) exp(-h^2))
#'
#' @export
simulate_spde <- function(nx = 101, nt = 101, sigma = 1,
                           lambda = 2, kappa = 0.5,
                           kernel = function(h) exp(-abs(h))) {
  # kernel should be a function of a single argument h = x - y

  f <- matrix(0, nx, nt)

  # Continuous part
  x_grid_aux <- seq(0, 2, length.out = 2 * nx)

  # Build the stationary covariance matrix
  q <- sigma * outer(x_grid_aux, x_grid_aux, Vectorize(function(x, y) kernel(x - y)))

  nx_aux <- length(x_grid_aux)

  noise <- MASS::mvrnorm(nt, mu = numeric(nx_aux), Sigma = q / nt)

  samples <- matrix(0, nt, nx_aux)
  for (i in 2:nt) {
    samples[i, 1:(nx_aux - i)] <- samples[i - 1, 2:(nx_aux + 1 - i)]
    samples[i, ] <- samples[i, ] + noise[i, ]
  }
  samples <- samples[1:nt, 1:nx]

  # Discontinuous (jump) part
  Int.Arr.times <- rexp(n = 100, rate = lambda)
  Arr.times <- cumsum(Int.Arr.times)
  Arr.times <- Arr.times[Arr.times < 1]

  jump.locations <- trunc(Arr.times * nt) + 1
  Jump.number <- length(Arr.times)

  cpp <- numeric(nt)
  if (Jump.number > 0) {
    if (1 %in% jump.locations) cpp[1] <- rnorm(1, sd = sqrt(kappa))
    for (i in 2:nt) {
      cpp[i] <- cpp[i - 1]
      if (i %in% jump.locations) cpp[i] <- cpp[i] + rnorm(1, sd = sqrt(kappa))
    }
  }

  CPP <- diag(cpp) %*% matrix(1, nt, nx)

  return(list(cont_samples = samples, discont_samples = samples + CPP, locs = jump.locations))
}
