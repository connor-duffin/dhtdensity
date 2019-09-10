#' Read the raw data into R
#'
#' @return A named list of the density measurements \code{y}, and the depths \code{z}.
read_data <- function() {
  dat_one <- Crux_KP150_Phs1
  dat_two <- Crux_KP150_Phs2
  dat_two$density_mat$`-11.7` <- NULL # remove these measurements
  dat_full <- rbind(dat_one$density_mat, dat_two$density_mat)
  y <- as.matrix(dat_full)
  dimnames(y) <- NULL
  z <- dat_one$depths

  list(y = y, z = z)
}
