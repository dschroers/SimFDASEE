#' Run Monte Carlo Simulation for BMTSM Functional Data
#'
#' This function performs a Monte Carlo simulation to evaluate different
#' covariance estimation methods for BMTSM (Brownian motion-type stochastic
#' models) functional data, optionally with jumps. It computes relative
#' errors of the estimated covariance matrices and effective dimensions
#' for multiple methods (raw, shifted/adjusted, and truncated versions).
#'
#' @param K Integer. Number of Monte Carlo iterations. Default is 10000.
#' @param nx Integer. Number of spatial/grid points. Default is 101.
#' @param nt Integer. Number of temporal points (observations). Default is 101.
#' @param kernel Function. Covariance kernel function taking a lag `h` as input. Default is exponential kernel: `function(h) exp(-abs(h))`.
#' @param sigma Numeric. Scale parameter for the simulation. Default is 1.
#' @param lambda Numeric. Jump intensity parameter for the simulation. Default is 2.
#' @param kappa Numeric. Jump size parameter for the simulation. Default is 0.5.
#' @param Delta Numeric. Grid spacing for truncation procedure. Default is 0.01.
#' @param sd_trunc Numeric. Standard deviation multiplier for truncation. Default is 3.
#' @param rho Numeric. Variance-explained threshold for effective dimension calculations. Default is 0.95.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{results}{A data frame containing relative errors and effective dimensions for different estimation methods: `rcv`, `sarcv`, `trunc_rcv`, `trunc_sarcv`. Each entry reports median with 0.25 and 0.75 quantiles in the format `median (Q1, Q3)`.}
#'   \item{d_true}{Effective dimension of the true covariance kernel based on `rho`.}
#'   \item{q}{The true covariance matrix used in the simulation.}
#' }
#'
#' @examples
#' # Run a small simulation with 10 iterations and a Gaussian kernel
#' gaussian_kernel <- function(h) exp(-h^2)
#' sim_res <- run_mc_simulation(K = 10, kernel = gaussian_kernel, nx = 50, nt = 51)
#' print(sim_res$results)
#'
#' @export
run_mc_simulation <- function(K = 10000,
                              nx = 101, nt = 101,
                              kernel = function(h) exp(-abs(h)),
                              sigma = 1, lambda = 2, kappa = 0.5,
                              Delta = 0.01, sd_trunc = 3,
                              rho = 0.95) {


  # Grid
  x_grid <- seq(0, 1, length.out = (nx-1))

  # True covariance matrix
  q <- outer(x_grid, x_grid, Vectorize(function(x, y) kernel(x - y)))
  E_q <- eigen(q)
  d_true <- which.max(cumsum(E_q$values) / sum(E_q$values) > rho)

  # Storage
  rel.err.rcv <- rel.err.sarcv <- rel.err.trunc_rcv <- rel.err.trunc_sarcv <- numeric(K)
  d_rv <- d_sarcv <- d_trunc_rv <- d_trunc_sarcv <- numeric(K)

  for (k in 1:K) {
    # Simulate data
    data <- simulate_spde(nx = nx, nt = nt, sigma = sigma, lambda = lambda, kappa = kappa, kernel = kernel)
    data_cont <- data$cont_samples
    data_discont <- data$discont_samples

    # Increments
    incr_cont <- data_cont[2:nt, 1:(nx-1)] - data_cont[1:(nt-1), 1:(nx-1)]
    adj_incr_cont <- data_cont[2:nt, 1:(nx-1)] - data_cont[1:(nt-1), 2:nx]

    incr_discont <- data_discont[2:nt, 1:(nx-1)] - data_discont[1:(nt-1), 1:(nx-1)]
    adj_incr_discont <- data_discont[2:nt, 1:(nx-1)] - data_discont[1:(nt-1), 2:nx]

    # Preliminary covariation
    C_prel_adj <- preliminary_covariation_estimator(data = adj_incr_discont, tq = 0.75)
    C_prel_reg <- preliminary_covariation_estimator(data = incr_discont, tq = 0.75)

    # Truncation locations
    locs_adj <- functional_data_truncation(d = d_star(C = C_prel_adj, rho = 0.75),
                                           C = C_prel_adj, data = adj_incr_discont, Delta = Delta, sd = sd_trunc)$locations
    locs_reg <- functional_data_truncation(d = d_star(C = C_prel_reg, rho = 0.75),
                                           C = C_prel_reg, data = incr_discont, Delta = Delta, sd = sd_trunc)$locations

    trunc_incr <- if (length(locs_reg) == 0) incr_discont else incr_discont[-locs_reg, ]
    trunc_adj_incr <- if (length(locs_adj) == 0) adj_incr_discont else adj_incr_discont[-locs_adj, ]

    # Covariances
    rv <- t(incr_cont) %*% incr_cont
    sarcv <- t(adj_incr_cont) %*% adj_incr_cont
    trunc_rv <- t(trunc_incr) %*% trunc_incr
    trunc_sarcv <- t(trunc_adj_incr) %*% trunc_adj_incr

    # Relative errors
    rel.err.rcv[k] <- hilbert_schmidt_norm(rv - q)
    rel.err.sarcv[k] <- hilbert_schmidt_norm(sarcv - q)
    rel.err.trunc_rcv[k] <- hilbert_schmidt_norm(trunc_rv - q)
    rel.err.trunc_sarcv[k] <- hilbert_schmidt_norm(trunc_sarcv - q)

    # Effective dimensions
    d_rv[k] <- which.max((cumsum(eigen(rv)$values) / sum(eigen(rv)$values)) > rho)
    d_sarcv[k] <- which.max((cumsum(eigen(sarcv)$values) / sum(eigen(sarcv)$values)) > rho)
    d_trunc_rv[k] <- which.max((cumsum(eigen(trunc_rv)$values) / sum(eigen(trunc_rv)$values)) > rho)
    d_trunc_sarcv[k] <- which.max((cumsum(eigen(trunc_sarcv)$values) / sum(eigen(trunc_sarcv)$values)) > rho)

    cat("Iterations left:", K - k, "\n")
  }

  # Compile results with median (Q1, Q3)
  rel.err.mat <- cbind(rel.err.rcv, rel.err.sarcv, rel.err.trunc_rcv, rel.err.trunc_sarcv) / hilbert_schmidt_norm(q)
  d.mat <- cbind(d_rv, d_sarcv, d_trunc_rv, d_trunc_sarcv)

  # Format relative errors
  rel_err_fmt <- sapply(1:4, function(i) {
    vals <- quantile(rel.err.mat[, i], probs = c(0.25, 0.5, 0.75))
    sprintf("%.3f (%.3f, %.3f)", vals[2], vals[1], vals[3])
  })

  # Format effective dimensions
  d_fmt <- sapply(1:4, function(i) {
    vals <- quantile(d.mat[, i], probs = c(0.25, 0.5, 0.75))
    sprintf("%.0f (%.0f, %.0f)", vals[2], vals[1], vals[3])
  })

  # Combine into data.frame
  results <- data.frame(
    rbind(rel_err = rel_err_fmt, d = d_fmt),
    stringsAsFactors = FALSE
  )
  colnames(results) <- c("rcv", "sarcv", "trunc_rcv", "trunc_sarcv")
  rownames(results) <- c("rel_err", "d")

  return(list(results = results, d_true = d_true, q = q))
}
