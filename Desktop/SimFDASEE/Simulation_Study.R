#library(SimFDASEE)

# Set seed for reproducibility
set.seed(123)

# Number of Monte Carlo iterations
K <- 10000

# Define scenarios
scenarios <- list(
  exponential = list(kernel = function(h) exp(-abs(h))),
  gaussian = list(kernel = function(h) exp(-h^2))
)

# Storage for results
simulation_results <- list()

# Run simulations for each scenario
for (scen_name in names(scenarios)) {
  cat("Running scenario:", scen_name, "\n")
  kernel_fun <- scenarios[[scen_name]]$kernel

  res <- run_mc_simulation(
    K = K,
    nx = 101,
    nt = 101,
    kernel = kernel_fun,
    sigma = 1,
    lambda = 2,
    kappa = 0.5,
    Delta = 0.01,
    sd_trunc = 3,
    rho = 0.95
  )

  simulation_results[[scen_name]] <- res
}


# Save each scenario's results
for (scen_name in names(simulation_results)) {
  file_path <- paste0("results/simulation_results_", scen_name, ".csv")
  write.csv(simulation_results[[scen_name]]$results,
            file = file_path,
            row.names = FALSE)
  cat("Saved results for scenario:", scen_name, "->", file_path, "\n")
}
