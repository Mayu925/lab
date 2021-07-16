source("common/gpcm_util_mayu.R")
source("common/uirt_util_mayu_5_25.R")

get_prm_list <- function(param, i, r, t){
  return(list(alpha_r = param$alpha_r[r], alpha_rt = param$alpha_rt, category_prm = append(0, param$beta_rk[r,])))
}

logit <- function(param, category_prm, x){
  return(param$alpha_r * ( x - param$alpha_rt - category_prm))
}

fisher_information <-function(param, theta) {
  Z <- fisher_information_no_alpha_multiplication(param, length(param$category_prm), theta)
  inf <- Z * param$alpha_i^2 * param$alpha_r^2  * param$alpha_c^2  * 1.7^2 
  return(inf)
}

get_param_size <- function(data){
  return(3* data$R + data$R * (data$K - 2) + data$T * 2  +data$I * 2 + data$J - (5 + data$T))
}

get_estimates <- function(fit1, setting){

  theta <- summary(fit1, par="theta")$summary[,"mean"]
  alpha_r <- summary(fit1, par="alpha_r")$summary[,"mean"]
  alpha_rt <- convert_alpha_rt(summary(fit1, par="alpha_rt")$summary[,"mean"] ,setting$n_time, setting$n_rater)
  category_prm <- summary(fit1, par="beta_rk")$summary[,"mean"]
  param = list(theta = theta, alpha_r = alpha_r, alpha_rt = alpha_rt, beta_rk = category_prm)
  return(param)
}

get_Rhat_stat <- function(fit1){
  RhatData <- c( 
                 summary(fit1, par="alpha_r")$summary[,"Rhat"],
                 summary(fit1, par="alpha_rt")$summary[,"Rhat"],
                 summary(fit1, par="beta_rk")$summary[,"Rhat"],
                 summary(fit1, par="theta")$summary[,"Rhat"])
  return(list(meanRhat = mean(RhatData), maxRhat = max(RhatData), countOver11 = sum(RhatData > 1.1)))
}

generate_true_param <- function(setting){
  theta <- rnorm(setting$n_person, 0, 1.0)
  alpha_r <- generate_constrained_alpha(setting$n_rater)
  alpha_rt <- rnorm(setting$n_rater*setting$n_time)
  beta_rk <- gen_category_param(setting$n_rater, setting$K)
  theta = theta - mean(theta)
  param = list(theta = theta, alpha_r = alpha_r, alpha_rt = alpha_rt,
               beta_rk = beta_rk)
  return(param)
}

get_error <- function(true_param, est_param){
  error_alpha_r <- calc_error(true_param$alpha_r, est_param$alpha_r)
  error_alpha_rt <- calc_error(true_param$alpha_rt, est_param$alpha_rt)
  error_category <- list(RMSE = sqrt(apply((true_param$beta_rk - est_param$beta_rk)^2, 2, mean)), 
                         BIAS = apply((true_param$beta_rk - est_param$beta_rk), 2, mean));
  error_theta <- calc_error(true_param$theta, est_param$theta)
  rmse <- list(theta = error_theta, alpha_r = error_alpha_r, alpha_rt = error_alpha_rt,
                beta_rk = error_category)  
  return(rmse)
}

