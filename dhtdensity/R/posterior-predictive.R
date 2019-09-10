#' Sample from the posterior predictive distribution
#'
#' @param mcmc_output A list output from \code{sample_dht_mcmc} or
#' \code{sample_dht_mcmc_ind}.
#' @param z_grid A grid of depths to evaluate the function over.
#' @param time A time to compute the posterior predictive distribution for.
#' @param n_rep Number of posterior predictive samples.
#' @param n_warmup Number of warmup iterations. Must not be on the global scale,
#' but on the actual number of total iterations scale.
#' @return A matrix of posterior predictive draws. Rows are individual density curves.
post_pred <- function(z_grid, n_rep, n_warmup, mcmc_output, time) {
  n_grid <- length(z_grid)
  y_rep <- matrix(NA, nrow = n_rep, ncol = n_grid)

  alpha <- mcmc_output$alpha
  beta <- mcmc_output$beta
  sigma2_y <- mcmc_output$sigma2_y

  for (i in 1:n_rep) {
    alpha_curr <- alpha[time, , n_warmup + i]
    beta_curr <- beta[time, , n_warmup + i]
    mean <- alpha_curr[1] + alpha_curr[2] * (
      tanh((z_grid + beta_curr[1]) / beta_curr[2]) + tanh((z_grid + beta_curr[3]) / beta_curr[4])
    )
    y_rep[i, ] <- rnorm(n_grid, mean, sigma2_y[n_warmup + i])
  }

  y_rep
}
