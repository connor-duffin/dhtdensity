dht <- function(z, beta) {
  tanh((z + beta[1]) / beta[2]) + tanh((z + beta[3]) / beta[4])
}

log_beta_post_ind <- function(beta, y, f, beta_prior_mean, beta_prior_sd,
                              alpha_curr, sigma2_curr) {
  sum(dnorm(y, mean = f %*% alpha_curr, sd = sqrt(sigma2_curr), log = TRUE)) +
    sum(dnorm(beta, beta_prior_mean, beta_prior_sd, log = TRUE))
}

#' MCMC sampler for non-hierarchical model (no hierarchy)
#'
#' @param y Vector of observed values
#' @param n_sample Number of MCMC samples.
#' @param n_adapt Number of adaptation samples
#' @param n_adapt_interval Interval between adaptations.
#' @param n_warmup Number of warmup iterations.
#' @param priors List of prior values. Can contain \code{beta_prior_mean} and
#' \code{beta_prior_sd}.
#' @param inits List of initial values. Can contain \code{alpha}, \code{beta},
#' and \code{sigma2_y}. If any unspecified reasonable random starting values 
#' will be generated.
#' @return A list of all MCMC samples of alpha, beta, and sigma2_y.
sample_dht_mcmc_ind <- function(
  y,
  z,
  n_sample,
  n_adapt,
  n_adapt_interval,
  n_warmup,
  chol_prop_mat,
  priors = list(
    alpha_mean = NULL,
    alpha_sd = NULL,
    beta_mean = NULL,
    beta_sd = NULL
  ),
  inits = list(
    alpha = NULL,
    beta = NULL,
    sigma2_y = NULL
  ),
  control = list(
    thin = 1,
    progress = FALSE
  )
) {
  stopifnot(!any(
    is.null(priors$alpha_mean),
    is.null(priors$alpha_sd),
    is.null(priors$beta_mean),
    is.null(priors$beta_sd)
  ))

  if (is.null(inits$alpha)) {
    inits$alpha <- c(1025, -5) + rnorm(2, sd = 0.05)
  }
  if (is.null(inits$beta)) {
    inits$beta <- c(75, 80, 150, 80) + rnorm(4, sd = 0.5)
  }
  if (is.null(inits$sigma2_y)) {
    inits$sigma2_y_curr <- rep(0.01, 4) + rnorm(4, sd = 0.001)
  }

  alpha_curr <- inits$alpha
  beta_curr <- inits$beta
  sigma2_y_curr <- inits$beta_curr

  n <- length(y)
  n_alpha <- length(alpha_curr)
  n_beta <- length(beta_curr)
  thin_index <- 1

  if ((n_sample - 1) %% control$thin == 0) {
    n_save <- (n_sample - 1) %/% control$thin + 1
  } else {
    n_save <- (n_sample - 1) %/% control$thin + 2
  }

  alpha <- matrix(NA, nrow = n_save, ncol = n_alpha)
  beta <- matrix(NA, nrow = n_save, ncol = n_beta)
  sigma2_y <- matrix(NA, nrow = n_save, ncol = 1)

  for (i in 1:n_sample) {
    if (control$progress) {
      cat("\rIteration: ", i, " / ", n_sample, sep = "")
    }
    f_curr <- cbind(1, dht(z, beta_curr))
    
    # sigma2_y update: gibbs
    sigma2_y_curr <- 1 / rgamma(
      1, n / 2, rate = sum((y - f_curr %*% alpha_curr)^2) / 2
    )
    # alpha update: gibbs
    alpha_cov <- solve(
      t(f_curr) %*% f_curr / sigma2_y_curr + diag(1 / priors$alpha_sd^2 + 1e-8)
    )
    alpha_mean <- alpha_cov %*%
      (diag(1 / priors$alpha_sd^2) %*% priors$alpha_mean + t(f_curr) %*% y / sigma2_y_curr)
    alpha_curr <- alpha_mean + chol(alpha_cov) %*% rnorm(n_alpha)

    if (i >= n_adapt & i < n_warmup & i %% n_adapt_interval == 0) {
      chol_prop_mat <- chol(
        2.38^2 * (cov(beta[1:(i - 1), ]) / 4 + 1e-8 * diag(n_beta))
      )
    }

    # beta update: mh step
    beta_prop <- beta_curr + chol_prop_mat %*% rnorm(n_beta) / 4
    f_prop <- cbind(1, dht(z, beta_prop))
    ratio_log <- min(
      log_beta_post_ind(
        beta_prop, y, f_prop, priors$beta_mean, priors$beta_sd, alpha_curr, sigma2_y_curr
      ) -
        log_beta_post_ind(
          beta_curr, y, f_curr, priors$beta_mean, priors$beta_sd, alpha_curr, sigma2_y_curr
        ),
      0
    ) 
    if (log(runif(1)) < ratio_log) {
      beta_curr <- beta_prop
    }

    # store sampled values
    if ((i - 1) %% control$thin == 0 || i == n_sample) {
      sigma2_y[thin_index, ] <- sigma2_y_curr
      alpha[thin_index, ] <- alpha_curr
      beta[thin_index, ] <- beta_curr
      thin_index <- thin_index + 1
    }
  }
  if (control$progress) {
    cat("\n")
  }

  alpha <- coda::mcmc(alpha)
  beta <- coda::mcmc(beta)
  sigma2_y <- coda::mcmc(sigma2_y)

  colnames(alpha) <- paste("alpha", 1:n_alpha, sep = "")
  colnames(beta) <- paste("beta", 1:n_beta, sep = "")
  colnames(sigma2_y) <- "sigma2_y"

  output <- list(
    alpha = alpha,
    beta = beta,
    sigma2_y = sigma2_y
  )

  output
}
