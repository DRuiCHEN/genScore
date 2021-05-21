library(foreach)

#' @param n_s,n_t sample size
#' @param d covariate dimension
#' @return a list containing x, a, y and oracle info
make_data <- function(n_s, n_t, d,
                      x_setting = 1,
                      ps_setting = 1,
                      rho_setting = 1,
                      y_setting = 1,
                      err_sigma = 1) {
  
  x_fn <- .make_x_fn(x_setting, d)
  ps_fn <- .make_prob_fn(ps_setting)
  rho_fn <- .make_prob_fn(rho_setting)
  y_fn <- .make_y_fn(y_setting)
  
  # covariates
  xs <- foreach(1:n_s, .combine = "rbind") %do% {
    repeat {
      x <- x_fn()
      accept <- rho_fn(x)
      if (runif(1) < accept) return(x)
    }
  }
  xt <- foreach(1:n_t, .combine = "rbind") %do% {
    repeat {
      x <- x_fn()
      accept <- 1 - rho_fn(x)
      if (runif(1) < accept) return(x)
    }
  }
  
  # true PS and treatment assignment
  p_true <- ps_fn(xs)
  a <- runif(n_s) <= p_true
  
  # observed outcome
  y_mat <- y_fn(xs)
  cate_true <- y_mat[, 2] - y_mat[, 1]
  y0_true <- y_mat[, 1]
  y <- y0_true + a * cate_true + rnorm(n_s) * err_sigma
  
  y_t_mat <- y_fn(xt)
  
  list(
    s = list(x = xs, a = a, y = y,
             oracle_info = data.frame(p_true, cate_true, y0_true)),
    t = list(x = xt,
             oracle_info = data.frame(p_true = ps_fn(xt), 
                                      cate_true = y_t_mat[, 2] - y_t_mat[, 1],
                                      y0_true = y_t_mat[, 1]))
  )
}




# Settings ----------------------------------------------------------------

.make_x_fn <- function(setting, d) {
  switch(
    setting,
    `1` = function(n = 1) matrix(rnorm(d), nrow = n)
  )
}

.make_prob_fn <- function(setting) {
  switch(
    setting,
    `1` = function(x) rep(.5, nrow(x)),
    `2` = function(x) plogis(.4 * x[,1] + .4 * x[,2] + .4 * x[,3]),
    `3` = function(x) plogis(.8 * x[,1] + .8 * x[,2] + .8 * x[,3]),
    `4` = function(x) plogis(.4 * x[,1] + .3 * x[,2]^3 + .2 * x[,3]^2),
    `5` = function(x) plogis(.8 * x[,1] + .6 * x[,2]^3 + .4 * x[,3]^2),
    `6` = function(x) plogis(.3 * x[,1] - .3 * x[,3]),
    `7` = function(x) plogis(.3 * x[,1] - .4 * x[,3]^2)
  ) 
}

.make_y_fn <- function(setting) {
  switch(
    setting,
    `1` = function(x) {
      y0 <- x[,1]
      y1 <- y0 +  .5 * x[,1] + x[,2] + 1
      cbind(y0, y1)
    },
    `2` = function(x) {
      mid <- x[,1] - x[,4]
      dif <- .4 * (x[,1] - .5)^2 + .5 * x[, 2]^2
      cbind(mid - .5 * dif, mid + .5 * dif)
    },
    `3` = function(x) {
      mid <- 1.2 + x[, 1] * .17 + x[, 2] * .16 - x[, 4] * .15
      dif <- (x[, 1] > -.5) * (x[, 2] < .5)
      y0 <- exp(mid - .3 * dif)
      y1 <- exp(mid - .3 * (1 - dif))
      cbind(y0, y1)
    }
  )
}
