get_prm_list <- function(param, i, r, t){
  return(list(
              beta_rt = param$beta_rt[r,t], 
              category_prm = append(0, param$beta_rk[r,])))
}

logit <- function(param, category_prm, x){
  return(x - param$beta_rt - category_prm)
}

get_estimates <- function(fit, setting){
  theta <- summary(fit, par="theta")$summary[,"mean"]
  beta_rt <- convert_beta_rt(summary(fit, par="beta_rt")$summary[,"mean"] ,setting$n_time, setting$n_rater)
  category_prm <- convert_category_estimates(summary(fit, par="beta_rk")$summary[,"mean"],setting$n_rater, setting$K)
  param = list(theta = theta, 
               beta_rt = beta_rt, 
               beta_rk = category_prm)
  return(param)
}

get_Rhat_stat <- function(fit){
  RhatData <- c( 
                 summary(fit, par="beta_rt")$summary[,"Rhat"],
                 summary(fit, par="beta_rk")$summary[,"Rhat"],
                 summary(fit, par="theta")$summary[,"Rhat"])
  return(list(meanRhat = mean(RhatData), maxRhat = max(RhatData), countOver11 = sum(RhatData > 1.1)))
}

generate_true_param <- function(setting){
  theta <- rnorm(setting$n_person, 0, 1.0)
  beta_rt <- generate_beta_rt_self(setting$n_rater, setting$n_time)
  beta_rk <- gen_category_param(setting$n_rater, setting$K)
  theta = theta - mean(theta)
  param = list(theta = theta, 
               beta_rt = beta_rt,
               beta_rk = beta_rk)
  return(param)
}

get_error <- function(true_param, est_param){
  error_beta_rt <- list(RMSE = sqrt(apply((true_param$beta_rt-est_param$beta_rt)^2, 2, mean)),
                         BIAS = apply((true_param$beta_rt - est_param$beta_rt), 2, mean));
  error_category <- list(RMSE = sqrt(apply((true_param$beta_rk - est_param$beta_rk)^2, 2, mean)), 
                         BIAS = apply((true_param$beta_rk - est_param$beta_rk), 2, mean));
  error_theta <- calc_error(true_param$theta, est_param$theta)
  rmse <- list(theta = error_theta, 
               beta_rt = error_beta_rt,
               beta_rk = error_category)
  return(rmse)
}
