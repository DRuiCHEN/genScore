#' Fit a linear or logistic regression and output predictions
#'
#' @param x covariates.
#' @param label labels, either numeric or logical.
#' @param x_extra additional covariates that we want to get predictions.
#' @return A list of predictions: pred(, pred_extra).
.fit <- function(x, label, x_extra = NULL)
{
  faml <- switch(typeof(label),
                 logical = "binomial",
                 double = "gaussian")
  
  mod <- suppressWarnings(glm(label ~ x, family = faml))
  
  .get_pred <- function(x) {
    coef <- mod$coefficients
    coef[is.na(coef)] <- 0
    z <- coef[1] + drop(x %*% coef[-1])
    if (is.logical(label)) z <- plogis(z)
    z
  }
  
  out <- list(pred = .get_pred(x))
  if (!is.null(x_extra)) out$pred_extra <- .get_pred(x_extra)
  out
}


#' Estimate the propensity score
#' 
#' @param xs covariates of the source sample.
#' @param as treatment assignment (binary) of the source sample.
#' @param xt covariates of the target sample.
#' @param x_test covariates of additional test data.
#' @return A list of estimated propensity scores: s(, t, test).
get_prop_score <- function(xs, as, xt = NULL, x_test = NULL)
{
  tmp <- .fit(xs, as == 1, rbind(xt, x_test))
  out <- list(s = tmp$pred)
  if (!is.null(xt)) out$t <- tmp$pred_extra[1:NROW(xt)]
  if (!is.null(x_test)) out$test <- tmp$pred_extra[NROW(xt) + 1:NROW(x_test)]
  out
}


#' Estimate the participation probability
#' 
#' @param xs,xt,x_test covariates.
#' @return A list of estimated participation probabilities: s, t(, test).
get_prob_source <- function(xs, xt, x_test = NULL)
{
  pred <- .fit(rbind(xs, xt), rep(c(TRUE, FALSE), c(NROW(xs), NROW(xt))), x_test)
  out <- list(s = pred$pred[1:NROW(xs)],
              t = pred$pred[NROW(xs) + 1:NROW(xt)])
  if (!is.null(x_test)) out$test <- pred$pred_extra
  out
}


#' Outcome regression
#'
#' @param xs,as,ys covariates, treatment assignment, observed outcome of the
#'   source sample.
#' @param xt covariates of the target sample.
#' @return A list of estimated potential outcomes: y1, y0, y1t, y0t.
get_mu <- function(xs, as, ys, xt)
{
  if (missing(xt)) {
    out <- list(y1 = .fit(xs[as == 1,], ys[as == 1], xs)$pred_extra,
                y0 = .fit(xs[as != 1,], ys[as != 1], xs)$pred_extra)
  }
  else{
    pred1 <- .fit(xs[as == 1,], ys[as == 1], rbind(xs, xt))$pred_extra
    pred0 <- .fit(xs[as != 1,], ys[as != 1], rbind(xs, xt))$pred_extra
    out <- list(y1 = pred1[1:NROW(xs)],
                y0 = pred0[1:NROW(xs)],
                y1t = pred1[NROW(xs) + 1:NROW(xt)],
                y0t = pred0[NROW(xs) + 1:NROW(xt)])
    
  }
  out
}


#' Estimate ATE using the IPW, outcome regression (OR) and augmented IPW (DR)
#' estimators
#'
#' @param x,a,y covariates, treatment assignment and observed outcome.
#' @param normalized_weight logical; whether to normalize the weights within the
#'   treated/control group.
#' @param se_boot logical; if TRUE, will perform bootstrap to obtain standard
#'   error estimates and estimation quantiles.
#' @param n_boot number of bootstrap resampling.
#' @param alpha the alpha/2 and 1-alpha/2 quantiles will be output.
#' @return IPW, OR, DR estimates. If se_boot is TRUE, then the standard error
#'   estimates and estimation quantiles will also be reported.
estimate_ATE <- function(x, a, y, 
                         normalized_weight = TRUE,
                         se_boot = FALSE, n_boot = 1000, alpha = .05)
{
  n <- NROW(x)
  
  ps <- get_prop_score(x, a)$s
  mu <- get_mu(x, a, y)
  
  # Compute the IPW, OR, DR estimates
  OR <- mean(mu$y1 - mu$y0)
  wt <- ifelse(a == 1, 1/ps, 1/(1-ps))
  if (normalized_weight) {
    wt[a == 1] <- wt[a == 1] / sum(wt[a == 1]) * n
    wt[a != 1] <- wt[a != 1] / sum(wt[a != 1]) * n
  }
  IPW <- mean(wt * (2 * (a == 1) - 1) * y)
  err <- y - ifelse(a == 1, mu$y1, mu$y0)
  DR <- mean(wt * (2 * (a == 1) - 1) * err) + OR
  
  res <- list(estimates = c(IPW = IPW, OR = OR, DR = DR))
  
  # bootstrap
  if (se_boot) {
    est_boot <- sapply(
      1:n_boot,
      function(b) {
        # We use the repeat loop to make sure the bootstrap sample contains
        # treated and control data
        repeat {
          idx <- sample(n, replace = TRUE)
          if (sum(a[idx]) > 1 && sum(a[idx]) < n-1) break
        }
        estimate_ATE(x[idx,], a[idx], y[idx], normalized_weight = normalized_weight)$estimates
      }
    )
    
    res$se_boot <- apply(est_boot, 1, sd)
    res$ql_boot <- apply(est_boot, 1, function(z) quantile(z, alpha/2))
    res$qu_boot <- apply(est_boot, 1, function(z) quantile(z, 1-alpha/2))
  }
  
  res
}


#' Estimate ATE for the target population using the IPW, outcome regression (OR)
#' and augmented IPW (DR) estimators
#'
#' @param xs,as,ys covariates, treatment assignment and observed outcome of the
#'   source sample.
#' @param xt covariates of the target sample.
#' @param normalized_weight logical; whether to normalize the weights within the
#'   treated/control group of the source sample.
#' @param se_boot logical; if TRUE, will perform bootstrap to obtain standard
#'   error estimates and estimation quantiles.
#' @param n_boot number of bootstrap resampling.
#' @param alpha the alpha/2 and 1-alpha/2 quantiles will be output.
#' @return IPW, OR, DR estimates. If se_boot is TRUE, then the standard error
#'   estimates and estimation quantiles will also be reported.
estimate_ATET <- function(xs, as, ys, xt,
                          normalized_weight = TRUE,
                          se_boot = FALSE, n_boot = 1000, alpha = .05)
{
  ns <- NROW(xs)
  nt <- NROW(xt)
  
  ps <- get_prop_score(xs, as)$s
  rho <- get_prob_source(xs, xt)$s
  mu <- get_mu(xs, as, ys, xt)
  
  # Compute the IPW, OR, DR estimates
  OR <- mean(mu$y1t - mu$y0t)
  wt <-   wt <- ifelse(as == 1, 1/ps, 1/(1-ps)) * (1-rho) / rho * ns / nt
  if (normalized_weight) {
    wt[as == 1] <- wt[as == 1] / sum(wt[as == 1]) * ns
    wt[as != 1] <- wt[as != 1] / sum(wt[as != 1]) * ns
  }
  IPW <- mean(wt * (2 * (as == 1) - 1) * ys)
  err <- ys - ifelse(as == 1, mu$y1, mu$y0)
  DR <- mean(wt * (2 * (as == 1) - 1) * err) + OR
  
  res <- list(estimates = c(IPW = IPW, OR = OR, DR = DR))
  
  
  if (se_boot) {
    est_boot <- sapply(
      1:n_boot,
      function(b) {
        # The resampling is on the combined sample (source + target). We use the
        # repeat loop to make sure the bootstrap sample contains treated,
        # control and target data
        repeat {
          idx <- sample(ns + nt, replace = TRUE)
          idxs <- idx[idx <= ns]
          idxt <- idx[idx > ns] - ns
          if (length(idxs) > 0 && length(idxt) > 0 && sum(as[idxs]) > 1 && sum(as[idxs]) < ns-1) break
        }
        estimate_ATET(xs[idxs,], as[idxs], ys[idxs], xt[idxt,],
                      normalized_weight = normalized_weight)$estimates
      }
    )
    
    res$se_boot <- apply(est_boot, 1, sd)
    res$ql_boot <- apply(est_boot, 1, function(z) quantile(z, alpha/2))
    res$qu_boot <- apply(est_boot, 1, function(z) quantile(z, 1-alpha/2))
  }
  
  res
}