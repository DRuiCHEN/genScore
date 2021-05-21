library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel()
source("simu_make_data.R")
source("utils.R")
source("method.R")

d <- 5
n_s <- 600; n_t <- 800; n_test <- 1e5
n_simu <- 3
n_boot <- 1000
rho_setting <- 2; ps_setting <- 6; y_setting <- 2

set.seed(222)

system.time({
  result <- foreach(1:n_simu) %dopar% {
    dat <- make_data(n_s, n_t, d, 
                     rho_setting = rho_setting, ps_setting = ps_setting, y_setting = y_setting)
    
    dat_test <- make_data(1, n_test, d, 
                          rho_setting = rho_setting, ps_setting = ps_setting, y_setting = y_setting)
    
    # full
    ATE_t_full <- mean(dat_test$t$oracle_info$cate_true)
    SATE_t_full <- mean(dat$t$oracle_info$cate_true)
    res_full <- estimate_ATET(dat$s$x, dat$s$a, dat$s$y, dat$t$x, se_boot = TRUE, n_boot = n_boot)
    
    evaluation <- tibble(
      ATE_t = ATE_t_full,
      SATE_t = SATE_t_full,
      estimate = res_full$estimates,
      se = res_full$se_boot,
      boot_l = res_full$ql_boot,
      boot_u = res_full$qu_boot,
      type = names(res_full$estimates),
      subset = FALSE
    )
    
    # subset
    ps <- get_prop_score(dat$s$x, dat$s$a, dat$t$x, dat_test$t$x)
    rho <- get_prob_source(dat$s$x, dat$t$x, dat_test$t$x)  
    sub_res <- get_subpop(rho, ps)
    sub_ind_test <- (1 - rho$test) / rho$test / ps$test / (1 - ps$test) <= sub_res$thres
    
    ATE_t_sub <- mean(dat_test$t$oracle_info$cate_true[sub_ind_test])
    SATE_t_sub <- mean(dat$t$oracle_info$cate_true[sub_res$t])
    
    res_sub <- estimate_ATET(dat$s$x[sub_res$s,], dat$s$a[sub_res$s], dat$s$y[sub_res$s], 
                             dat$t$x[sub_res$t,], se_boot = TRUE, n_boot = n_boot)
    
    evaluation <- 
      bind_rows(
        evaluation,
        tibble(
          ATE_t = ATE_t_sub,
          SATE_t = SATE_t_sub,
          estimate = res_sub$estimates,
          se = res_sub$se_boot,
          boot_l = res_sub$ql_boot,
          boot_u = res_sub$qu_boot,
          type = names(res_sub$estimates),
          subset = TRUE
        )
      )
    
    list(evaluation = evaluation,
         sub_summary = c(prop_s = mean(sub_res$s),
                         prop_t = mean(sub_res$t),
                         prop_test = mean(sub_ind_test),
                         thres = sub_res$thres))
  }
})
