source("common/gpcm_util_mayu.R")
source("common/uirt_util_mayu_5_25.R")

get_prm_list <- function(param, i, r, t){
  return(list(alpha_r = param$alpha_r[r], pai_1r = param$pai_1r[r], pai_0r = param$pai_0r[r], category_prm = append(0, param$beta_rk[r,])))
}

logit <- function(param, category_prm, x){
  return(param$alpha_r * ( x - param$pai_0r - param$pai_1r * param$alpha_rt - category_prm))
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
  #alpha_i <- get_alpha_estimates_with_restriction(summary(fit, par="alpha_i")$summary[,"mean"])
  #alpha_c <- get_alpha_estimates_with_restriction(summary(fit, par="alpha_c")$summary[,"mean"])
  #beta_i <- get_estimates_with_mean_restriction(summary(fit, par="beta_i")$summary[,"mean"])
  #beta_r <- summary(fit, par="beta_r")$summary[,"mean"]
  #beta_c <- get_estimates_with_mean_restriction(summary(fit, par="beta_c")$summary[,"mean"])
  #tau_r <- get_alpha_estimates_with_restriction(summary(fit, par="tau_r")$summary[,"mean"])
  
  theta <- summary(fit1, par="theta")$summary[,"mean"]
  alpha_r <- summary(fit1, par="alpha_r")$summary[,"mean"]
  alpha_rt <- get_estimates_with_mean_restriction(summary(fit1, par="alpha_rt")$summary[,"mean"])
  category_prm <- convert_category_estimates(summary(fit1, par="beta_rk")$summary[,"mean"], setting$n_time, setting$K)
  pai_0r <- get_estimates_with_mean_restriction(summary(fit1, par="pai_0r")$summary[,"mean"])
  pai_1r <- get_alpha_estimates_with_restriction(summary(fit1, par="pai_1r")$summary[,"mean"])
  
  param = list(theta = theta, alpha_r = alpha_r, alpha_rt = alpha_rt, beta_rk = category_prm, pai_0r = pai_0r, pai_1r = pai_1r)
  return(param)
}

get_Rhat_stat <- function(fit1){
  RhatData <- c( 
                 summary(fit1, par="alpha_r")$summary[,"Rhat"],
                 summary(fit1, par="alpha_rt")$summary[,"Rhat"],
                 summary(fit1, par="beta_rk")$summary[,"Rhat"],
                 summary(fit1, par="pai_0r")$summary[,"Rhat"],
                 summary(fit1, par="pai_1r")$summary[,"Rhat"],
                 summary(fit1, par="theta")$summary[,"Rhat"])
  return(list(meanRhat = mean(RhatData), maxRhat = max(RhatData), countOver11 = sum(RhatData > 1.1)))
}

generate_true_param <- function(setting){
  theta <- rnorm(setting$n_person, 0, 1.0)
  #alpha_i <- generate_constrained_alpha(setting$n_item) 
  alpha_r <- generate_constrained_alpha(setting$n_rater)
  #alpha_c <- generate_constrained_alpha(setting$n_rubric)
  #beta_i <- rnorm(setting$n_item, 0, 1.0)
  #beta_r <- rnorm(setting$n_rater, 0, 1.0)
  #beta_c <- rnorm(setting$n_rubric, 0, 0.1)
  #tau_r <- generate_constrained_alpha(setting$n_rater) * 0.5 + 0.5
  alpha_rt <-
  beta_rk <- gen_category_param(setting$n_rater, setting$K)
  pai_0r <-
  pai_1r <-
  theta = theta - mean(theta)
  #beta_i = beta_i - mean(beta_i)
  #beta_r = beta_r - mean(beta_r)
  #beta_c = beta_c - mean(beta_c)
  param = list(theta = theta, alpha_r = alpha_r, alpha_rt = alpha_rt,
               beta_rk = beta_rk, pai_0r = pai_0r, pai_1r = pai_1r)
  return(param)
}

get_error <- function(true_param, est_param){
  error_alpha_i <- calc_error(true_param$alpha_i, est_param$alpha_i)
  error_alpha_r <- calc_error(true_param$alpha_r, est_param$alpha_r)
  error_alpha_c <- calc_error(true_param$alpha_c, est_param$alpha_c)
  error_beta_i <- calc_error(true_param$beta_i, est_param$beta_i)
  error_beta_r <- calc_error(true_param$beta_r, est_param$beta_r)
  error_beta_c <- calc_error(true_param$beta_c, est_param$beta_c)
  error_tau_r <- calc_error(true_param$tau_r, est_param$tau_r)
  error_category <- list(RMSE = sqrt(apply((true_param$beta_ck - est_param$beta_ck)^2, 2, mean)), 
                         BIAS = apply((true_param$beta_ck - est_param$beta_ck), 2, mean));
  error_theta <- calc_error(true_param$theta, est_param$theta)
  rmse <- list(theta = error_theta, alpha_i = error_alpha_i, alpha_r = error_alpha_r, alpha_c = error_alpha_c,
               beta_i = error_beta_i, beta_r = error_beta_r, beta_c = error_beta_c, 
               tau_r = error_tau_r, beta_ck = error_category)  
  return(rmse)
}

