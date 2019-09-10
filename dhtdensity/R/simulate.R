#' Simulate data according to the dht model
#'
#' @param z Depths at which to evaluate the density.
#' @param n_groups Number of groups to simulate.
#' @param sigma2_y Variance about the mean curve.
#' @param mu_alpha Mean value of the alpha's.
#' @param sigma2_alpha Variance of the alpha's.
#' @param mu_beta Mean value of the beta's.
#' @param sigma2_beta Variance of the beta's.
#' @return A matrix of simulated data. Rows are the individual density curves.
simulate_double_tanh <- function(z, n_groups, sigma2_y, mu_alpha, sigma2_alpha,
                                 mu_beta, sigma2_beta) {
  y <- matrix(NA, nrow = n_groups, ncol = length(z))

  for (i in seq(n_groups)) {
    alpha <- rnorm(length(mu_alpha), mu_alpha, sqrt(sigma2_alpha))
    beta <- rnorm(length(mu_beta), mu_beta, sqrt(sigma2_beta))
    y[i, ] <- rnorm(length(z),
      alpha[1] + alpha[2] * (
        tanh((z + beta[1]) / beta[2]) +
        tanh((z + beta[3]) / beta[4])
      ),
      sqrt(sigma2_y)
    )
  }

  y
}
