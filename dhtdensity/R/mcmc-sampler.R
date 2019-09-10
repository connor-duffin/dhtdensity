dht <- function(z, beta) {
  tanh((z + beta[1]) / beta[2]) + tanh((z + beta[3]) / beta[4])
}

log_beta_post <- function(beta, y, f, alpha_curr, sigma2_y_curr, mu_beta_curr,
                          sigma2_beta_curr) {
  sum(dnorm(y, mean = f %*% alpha_curr, sd = sqrt(sigma2_y_curr), log = TRUE)) +
    sum(dnorm(beta, mu_beta_curr, sqrt(sigma2_beta_curr), log = TRUE))
}

update_alpha <- function(n_alpha, y, f, sigma2_y, mu_alpha, sigma2_alpha) {
  alpha_cov <- solve(
    t(f) %*% f / sigma2_y + diag(1 / sigma2_alpha + 1e-8)
  )
  alpha_mean <- alpha_cov %*% 
    (diag(1 / sigma2_alpha) %*% mu_alpha + t(f) %*% y / sigma2_y)

  alpha_mean + chol(alpha_cov) %*% rnorm(n_alpha)
}

update_mu_alpha <- function(n_groups, n_alpha, alpha, sigma2_alpha, 
                            mu_alpha_prior_mean, mu_alpha_prior_sd) {
  mean <- (
    (colSums(alpha) / sigma2_alpha + mu_alpha_prior_mean / mu_alpha_prior_sd^2) /
    (n_groups / sigma2_alpha + 1 / mu_alpha_prior_sd^2)
  )
  sd <- sqrt(1 / (n_groups / sigma2_alpha + 1 / mu_alpha_prior_sd^2))

  rnorm(n_alpha, mean, sd) 
}

update_mu_beta <- function(n_groups, n_beta, beta, sigma2_beta,
                           mu_beta_prior_mean, mu_beta_prior_sd) {
  mean <- (
    (colSums(beta) / sigma2_beta + mu_beta_prior_mean / mu_beta_prior_sd^2) /
    (n_groups / sigma2_beta + 1 / mu_beta_prior_sd^2)
  )
  sd <- sqrt(1 / (n_groups / sigma2_beta + 1 / mu_beta_prior_sd^2))

  rnorm(n_beta, mean, sd) 
}

#' MCMC sampler for the DHT hierarchical model
#'
#' @param y Density data. Must be a matrix with rows corresponding to individual
#' groups.
#' @param z Grid of depths. Must be a vector of length the same as \code{ncol(y)}.
#' @param n_sample Number of MCMC samples.
#' @param chol_prop_mat An array of cholesky's of the desired proposal matrix for the
#' MH step for the beta's (nonlinear dht parameters).
#' @param inits A list of initial values. Can contain \code{alpha}, \code{beta},
#' \code{mu_alpha}, \code{mu_beta}, \code{sigma2_alpha}, \code{sigma2_beta}, 
#' and \code{sigma2_y}. Any values which are unspecified will take the
#' default values.
#' @param priors A list of prior values. Can contain \code{mu_alpha_mean}, \code{mu_alpha_mean},
#' \code{mu_beta_sd}, \code{mu_beta_sd}, and \code{sigma2_max}. Any values which
#' are unspecified will take the defaults.
#' @param control A list of control values. Can contain \code{thin}, and \code{progress}.
#' @return A named list containing elements \code{alpha}, \code{beta},
#' \code{mu_alpha}, \code{mu_beta}, \code{sigma2_alpha}, \code{sigma2_beta},
#' and \code{sigma2_y}.
sample_dht_mcmc <- function(
  y,
  z,
  n_sample,
  chol_prop_mat,
  inits = list(
    alpha = NULL,
    beta = NULL,
    mu_alpha = NULL,
    mu_beta = NULL,
    sigma2_alpha = NULL,
    sigma2_beta = NULL,
    sigma2_y = NULL
  ),
  priors = list(
    mu_alpha_mean = NULL,
    mu_beta_mean = NULL,
    mu_alpha_sd = NULL,
    mu_beta_sd = NULL,
    sigma2_max = 5000
  ),
  control = list(
    thin = 1,
    progress = FALSE
  )
) {
  n_groups <- nrow(y)
  n_y <- ncol(y)

  # these values crucial: must be specified
  stopifnot(!any(
    is.null(inits$alpha),
    is.null(inits$beta),
    is.null(inits$mu_alpha),
    is.null(inits$mu_beta)
  ))

  # these values are less finicky
  if (is.null(inits$sigma2_alpha)) {
    inits$sigma2_alpha <- rnorm(2, 1, 0.1)
  }
  if (is.null(inits$sigma2_beta)) {
    inits$sigma2_beta <- rnorm(4, 2, 0.1)
  }
  if (is.null(inits$sigma2_y)) {
    inits$sigma2_y <- rnorm(1, 0.01, 0.001)
  }

  alpha_curr <- inits$alpha
  beta_curr <- inits$beta
  mu_alpha_curr <- inits$mu_alpha
  mu_beta_curr <- inits$mu_beta
  sigma2_alpha_curr <- inits$sigma2_alpha
  sigma2_beta_curr <- inits$sigma2_beta
  sigma2_y_curr <- inits$sigma2_y

  n_alpha <- ncol(alpha_curr)
  n_beta <- ncol(beta_curr)
  thin_index <- 1

  if (is.matrix(chol_prop_mat)) {
    chol_prop_mat <- array(chol_prop_mat, dim = c(n_beta, n_beta, n_groups))
  }

  # these values MUST be specified
  stopifnot(!any(
    is.null(priors$mu_alpha_mean),
    is.null(priors$mu_beta_mean),
    is.null(priors$mu_alpha_sd),
    is.null(priors$mu_beta_sd)
  ))

  if (is.null(priors$sigma2_max)) {
    priors$sigma2_max <- 5000
  }

  if ((n_sample - 1) %% control$thin == 0) {
    n_save <- (n_sample - 1) %/% control$thin + 1
  } else {
    n_save <- (n_sample - 1) %/% control$thin + 2
  }

  alpha <- array(NA, dim = c(n_groups, n_alpha, n_save))
  beta <- array(NA, dim = c(n_groups, n_beta, n_save))
  mu_alpha <- matrix(NA, nrow = n_save, ncol = n_alpha)
  mu_beta <- matrix(NA, nrow = n_save, ncol = n_beta)
  sigma2_alpha <- matrix(NA, nrow = n_save, ncol = n_alpha)
  sigma2_beta <- matrix(NA, nrow = n_save, ncol = n_beta)
  sigma2_y <- matrix(NA, nrow = n_save, ncol = 1)

  colnames(alpha) <- paste("alpha", 1:n_alpha, sep = "_")
  colnames(beta) <- paste("beta", 1:n_beta, sep = "_")
  colnames(mu_alpha) <- paste("mu_alpha", 1:n_alpha, sep = "_")
  colnames(mu_beta) <- paste("mu_beta", 1:n_beta, sep = "_")
  colnames(sigma2_alpha) <- paste("sigma2_alpha", 1:n_alpha, sep = "_")
  colnames(sigma2_beta) <- paste("sigma2_beta", 1:n_beta, sep = "_")
  colnames(sigma2_y) <- "sigma2_y"

  f_curr <- cbind(1, dht(z, beta_curr[1, ]))
  f_prop <- cbind(1, dht(z, beta_curr[1, ]))

  for (i in 1:n_sample) {
    if (control$progress) {
      cat("\rIteration: ", i, " / ", n_sample, sep = "")
    }
    rate_sigma2_alpha <- 0
    rate_sigma2_beta <- 0
    rate_sigma2_y <- 0

    for (j in 1:n_groups) {
      f_curr[, 2] <- dht(z, beta_curr[j, ])

      # update the alpha's: gibbs step
      alpha_curr[j, ] <- update_alpha(
        n_alpha, y[j, ], f_curr, sigma2_y_curr, mu_alpha_curr, sigma2_alpha_curr
      )

      # update the beta's: mh step
      beta_prop <- beta_curr[j, ] + chol_prop_mat[, , j] %*% rnorm(n_beta) / 4
      f_prop[, 2] <- dht(z, beta_prop)
      ratio_log <- min(
        log_beta_post(
          beta_prop, y[j, ], f_prop, alpha_curr[j, ], sigma2_y_curr,
          mu_beta_curr, sigma2_beta_curr
        ) -
          log_beta_post(
            beta_curr[j, ], y[j, ], f_curr, alpha_curr[j, ], sigma2_y_curr,
            mu_beta_curr, sigma2_beta_curr
          ),
        0
      )
      if (log(runif(1)) < ratio_log) {
        beta_curr[j, ] <- beta_prop
        f_curr <- f_prop
      }

      # in case of label switching
      if (beta_curr[j, 1] > beta_curr[j, 3]) {
        temp <- beta_curr[j, 1]
        beta_curr[j, 1] <- beta_curr[j, 3]
        beta_curr[j, 3] <- temp
      }

      # update the rates for each parameter
      rate_sigma2_alpha <- rate_sigma2_alpha +
        (alpha_curr[j, ] - mu_alpha_curr)^2 / 2
      rate_sigma2_beta <- rate_sigma2_beta +
        (beta_curr[j, ] - mu_beta_curr)^2 / 2
      rate_sigma2_y <- rate_sigma2_y +
        sum((y[j, ]  - f_curr %*% alpha_curr[j, ])^2) / 2
    }

    # update each sigma^2: gibbs steps
    sigma2_alpha_prop <- 1 / rgamma(
      n_alpha, shape = n_groups / 2 - 1, rate = rate_sigma2_alpha
    )
    if (any(sigma2_alpha_prop <= priors$sigma2_max)) {
      sigma2_alpha_curr <- sigma2_alpha_prop
    }

    sigma2_beta_prop <- 1 / rgamma(
      n_beta, shape = n_groups / 2  - 1, rate = rate_sigma2_beta
    )
    if (any(sigma2_beta_prop <= priors$sigma2_max)) {
      sigma2_beta_curr <- sigma2_beta_prop
    }

    sigma2_y_prop <- 1 / rgamma(
      1, shape = n_y * n_groups / 2, rate = rate_sigma2_y
    ) 
    if (any(sigma2_y_prop <= priors$sigma2_max)) {
      sigma2_y_curr <- sigma2_y_prop
    }

    mu_alpha_curr <- update_mu_alpha(
      n_groups, n_alpha, alpha_curr, sigma2_alpha_curr, priors$mu_alpha_mean, priors$mu_alpha_sd
    )
    mu_beta_curr <- update_mu_beta(
      n_groups, n_beta, beta_curr, sigma2_beta_curr, priors$mu_beta_mean, priors$mu_beta_sd
    )

    # save the samples
    if ((i - 1) %% control$thin == 0 || i == n_sample) {
      alpha[, , thin_index] <- alpha_curr
      beta[, , thin_index] <- beta_curr

      mu_alpha[thin_index, ] <- mu_alpha_curr
      mu_beta[thin_index, ] <- mu_beta_curr

      sigma2_alpha[thin_index, ] <- sigma2_alpha_curr
      sigma2_beta[thin_index, ] <- sigma2_beta_curr
      sigma2_y[thin_index] <- sigma2_y_curr

      thin_index <- thin_index + 1
    }
  }
  if (control$progress) {
    cat("\n")
  }

  output <- list(
    alpha = acoda::mcmca(alpha, thin = control$thin),
    beta = acoda::mcmca(beta, thin = control$thin),
    mu_alpha = coda::mcmc(mu_alpha, thin = control$thin),
    mu_beta = coda::mcmc(mu_beta, thin = control$thin),
    sigma2_alpha = coda::mcmc(sigma2_alpha, thin = control$thin),
    sigma2_beta = coda::mcmc(sigma2_beta, thin = control$thin),
    sigma2_y = coda::mcmc(sigma2_y, thin = control$thin)
  )

  output
}
