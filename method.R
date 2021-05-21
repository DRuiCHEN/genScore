#' Get optimal subset
#'
#' @param rho list of participation probabilities, output of the get_source_prob
#'   function.
#' @param ps list of propensity scores, output of the get_prop_score function.
#' @param y_sigma_sq a numeric vector representing the noise variance of each
#'   units. Default as 1.
#' @param thres threshold value for trimming If not supplied, will be set as the
#'   one such that the maximum kappa in the subset is no greater than twice the
#'   average value.
#' @param max_round maximum number of iterations when searching the trimming
#'   threshold.
#' @param eps precision tolerance when search the trimming threshold.
#' @param more_target proportion of additional target population to include.
#' @param path a logical indicating whether to output a series of thresholds and
#'   the corresponding V(B).
#' @param echo_iter a logical indicating whether to echo during the iteration
#'   process.
#' @return A list of s: a logical vector indicating whether to include the
#'   source units; t: a logical vector indicating whether to include the target
#'   units; thres: the threshold used in determining the subpopulation; thr: a
#'   vector of thresholds; avar: a vector of V(B) estimates corresponding to
#'   thr.
get_subpop <- function(rho, ps, 
                       y_sigma_sq = 1,
                       thres, 
                       max_round = 2000, eps = 1e-8,
                       more_target = 0,
                       path = FALSE, echo_iter = FALSE) 
{
  
  rho_all <- c(rho$s, rho$t)
  ps_all <- c(ps$s, ps$t)
  n_s <- length(rho$s)
  n_t <- length(rho$t)
  
  w <- 1 - rho_all
  kappa_all <- (1 - rho_all) / rho_all / ps_all / (1 - ps_all) * y_sigma_sq
  
  # If the threshold is not supplied, find one using the following iteration:
  # thres_new = E[kappa | kappa <= thres_old, S = 0]
  if (missing(thres)) {
    thres <- Inf
    for (i in 1:max_round) {
      if (echo_iter) cat(sprintf("iteration %d: thres = %.3f\n", i, thres))
      old_thres <- thres
      sub_idx <- kappa_all <= thres
      thres <- 2 * weighted.mean(kappa_all[sub_idx], w[sub_idx])
      if (abs(thres - old_thres) < eps) break
      if (i == max_round) warning("Max iteration reached without convergence!")
    }
  }
  
  # If more_target is greater than 0, use a larger threshold
  if (more_target > 0) {
    thr <- sort(kappa_all[n_s + 1:n_t])
    n_tmp <- sum(thr <= thres)
    thres <- thr[floor(n_tmp + more_target * (n_t - n_tmp))]
  }
  
  # Subjects with kappa <= threshold are selected into the subset
  sub_idx <- kappa_all <= thres
  out <- list(s = sub_idx[1:n_s],
              t = sub_idx[n_s + 1:n_t],
              thres = thres)
  
  # Calculate asymptotic var vs different threshold values
  if (path) {
    out$thr <- sort(kappa_all[n_s + 1:n_t])
    out$avar <- sapply(out$thr, function(x) {
      ind <- kappa_all <= x
      mean((1 - rho_all) * kappa_all * ind) / mean((1 - rho_all) * ind)^2
    })
  }
  
  out
}
