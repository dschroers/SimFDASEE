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
#' @param rank_method True if the preliminary covariance estimator should be picked by a rank method.
#' @param w If preliminary truncation is conservative, the threshold is Delta^{w}.
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
                              rho = 0.95,
                              Trunc_rule = c("simple", "simple", "rank"),
                              w = c(.2, .4, .49)) {

  # Basic checks
  U <- length(Trunc_rule)
  if (length(w) != U) stop("`w` must have the same length as `Trunc_rule`.")

  # Grid (matches your increment dimensions: nx-1)
  x_grid <- seq(0, 1, length.out = (nx - 1))

  # True covariance matrix
  q <- outer(x_grid, x_grid, Vectorize(function(x, y) kernel(x - y)))
  E_q <- eigen(q)
  d_true <- which.max(cumsum(E_q$values) / sum(E_q$values) > rho)

  # Storage
  rel.err.rcv <- rel.err.sarcv <- numeric(K)
  d_rv <- d_sarcv <- numeric(K)

  rel.err.trunc_rcv   <- matrix(NA_real_, nrow = K, ncol = U)
  rel.err.trunc_sarcv <- matrix(NA_real_, nrow = K, ncol = U)
  d_trunc_rv          <- matrix(NA_real_, nrow = K, ncol = U)
  d_trunc_sarcv       <- matrix(NA_real_, nrow = K, ncol = U)

  for (k in seq_len(K)) {

    # Simulate data
    data <- simulate_spde(nx = nx, nt = nt, sigma = sigma, lambda = lambda, kappa = kappa, kernel = kernel)
    data_cont <- data$cont_samples
    data_discont <- data$discont_samples

    # Increments (continuous)
    incr_cont <- data_cont[2:nt, 1:(nx-1)] - data_cont[1:(nt-1), 1:(nx-1)]
    adj_incr_cont <- data_cont[2:nt, 1:(nx-1)] - data_cont[1:(nt-1), 2:nx]

    # Increments (discontinuous)
    incr_discont <- data_discont[2:nt, 1:(nx-1)] - data_discont[1:(nt-1), 1:(nx-1)]
    adj_incr_discont <- data_discont[2:nt, 1:(nx-1)] - data_discont[1:(nt-1), 2:nx]

    # Covariances (raw)
    rv <- t(incr_cont) %*% incr_cont
    sarcv <- t(adj_incr_cont) %*% adj_incr_cont

    # Relative errors (unscaled; scale later)
    rel.err.rcv[k] <- hilbert_schmidt_norm(rv - q)
    rel.err.sarcv[k] <- hilbert_schmidt_norm(sarcv - q)

    # Effective dimensions
    ev_rv <- eigen(rv, symmetric = TRUE, only.values = TRUE)$values
    ev_sarcv <- eigen(sarcv, symmetric = TRUE, only.values = TRUE)$values
    d_rv[k] <- which.max((cumsum(ev_rv) / sum(ev_rv)) > rho)
    d_sarcv[k] <- which.max((cumsum(ev_sarcv) / sum(ev_sarcv)) > rho)

    for (u in seq_len(U)) {

      if (Trunc_rule[u] == "rank") {

        # Preliminary covariation
        C_prel_adj <- preliminary_covariation_estimator(data = adj_incr_discont, tq = 0.75)
        C_prel_reg <- preliminary_covariation_estimator(data = incr_discont, tq = 0.75)

        # Truncation locations
        locs_adj <- functional_data_truncation(
          d = d_star(C = C_prel_adj, rho = 0.75),
          C = C_prel_adj, data = adj_incr_discont, Delta = Delta, sd = sd_trunc
        )$locations

        locs_reg <- functional_data_truncation(
          d = d_star(C = C_prel_reg, rho = 0.75),
          C = C_prel_reg, data = incr_discont, Delta = Delta, sd = sd_trunc
        )$locations

        trunc_incr <- if (length(locs_reg) == 0) incr_discont else incr_discont[-locs_reg, , drop = FALSE]
        trunc_adj_incr <- if (length(locs_adj) == 0) adj_incr_discont else adj_incr_discont[-locs_adj, , drop = FALSE]

        trunc_rv <- t(trunc_incr) %*% trunc_incr
        trunc_sarcv <- t(trunc_adj_incr) %*% trunc_adj_incr

      } else {
        # "simple" conservative truncation with threshold Delta^(2*w[u])
        norm_squares_adj <- rowSums(adj_incr_discont^2) * Delta
        norm_squares_reg <- rowSums(incr_discont^2) * Delta

        thr <- Delta^(2 * w[u])

        data_trunc_adj <- adj_incr_discont[norm_squares_adj < thr, , drop = FALSE]
        data_trunc_reg <- incr_discont[norm_squares_reg < thr, , drop = FALSE]

        trunc_sarcv <- t(data_trunc_adj) %*% data_trunc_adj
        trunc_rv <- t(data_trunc_reg) %*% data_trunc_reg
      }

      # Relative errors
      rel.err.trunc_rcv[k, u] <- hilbert_schmidt_norm(trunc_rv - q)
      rel.err.trunc_sarcv[k, u] <- hilbert_schmidt_norm(trunc_sarcv - q)

      # Effective dimensions
      ev_trv <- eigen(trunc_rv, symmetric = TRUE, only.values = TRUE)$values
      ev_tsarcv <- eigen(trunc_sarcv, symmetric = TRUE, only.values = TRUE)$values

      d_trunc_rv[k, u] <- which.max((cumsum(ev_trv) / sum(ev_trv)) > rho)
      d_trunc_sarcv[k, u] <- which.max((cumsum(ev_tsarcv) / sum(ev_tsarcv)) > rho)
    }

    cat("Iterations left:", K - k, "\n")
  }

  # ---------- Polishing / results table (NO V1/V2; keeps special chars) -------

  # Desired column names
  col_trunc_rcv   <- paste0("rcv_-(",   seq_len(U), ")")
  col_trunc_sarcv <- paste0("sarcv_-(", seq_len(U), ")")

  # Force colnames on matrices so they never become V1, V2, ...
  colnames(rel.err.trunc_rcv)   <- col_trunc_rcv
  colnames(rel.err.trunc_sarcv) <- col_trunc_sarcv
  colnames(d_trunc_rv)          <- col_trunc_rcv
  colnames(d_trunc_sarcv)       <- col_trunc_sarcv

  # Assemble matrices
  rel.err.mat <- cbind(
    rcv   = rel.err.rcv,
    sarcv = rel.err.sarcv,
    rel.err.trunc_rcv,
    rel.err.trunc_sarcv
  ) / hilbert_schmidt_norm(q)

  d.mat <- cbind(
    rcv   = d_rv,
    sarcv = d_sarcv,
    d_trunc_rv,
    d_trunc_sarcv
  )

  # Format helpers
  fmt_q <- function(x, digits = 3) {
    qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
            qs[2], qs[1], qs[3])
  }
  fmt_q_int <- function(x) {
    qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    sprintf("%.0f (%.0f, %.0f)", qs[2], qs[1], qs[3])
  }

  # Preserve names exactly (including '-' and parentheses)
  rel_err_fmt <- vapply(as.data.frame(rel.err.mat, check.names = FALSE), fmt_q, character(1))
  d_fmt       <- vapply(as.data.frame(d.mat,       check.names = FALSE), fmt_q_int, character(1))

  results <- as.data.frame(
    rbind(rel_err = rel_err_fmt, d = d_fmt),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  return(list(results = results, d_true = d_true, q = q))
}
