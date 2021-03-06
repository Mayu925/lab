get_prm_list <- function(param, r, t){
  return(list(pai_1r = param$pai_1r[r], 
    pai_0r = param$pai_0r[r], category_prm = append(0, param$beta_rk)))
}

logit <- function(setting, param, category_prm, x, t){
  return(x - param$pai_0r -param$pai_1r *(t /setting$n_time)- category_prm)
}

get_param_size <- function(data){
  return(3* data$R + data$R * (data$K - 2) + data$T * 2  +data$I * 2 + data$J - (5 + data$T))
}

get_estimates <- function(fit, setting){
  theta <- summary(fit, par="theta")$summary[,"mean"]
  #beta_rt <- convert_beta_rt(summary(fit, par="beta_rt")$summary[,"mean"] , setting$n_time, setting$n_rater)
  category_prm <- convert_category_estimates_previous(summary(fit, par="beta_rk")$summary[,"mean"],setting$K)
  pai_0r <- convert_pai1r(summary(fit, par="pai_0r")$summary[,"mean"], setting$n_rater)
  pai_1r <- summary(fit, par="pai_1r")$summary[,"mean"]
  param = list(theta = theta, 
               #beta_rt = beta_rt, 
               beta_rk = category_prm, 
               pai_1r = pai_1r,
               pai_0r = pai_0r
               )
  return(param)
}

get_Rhat_stat <- function(fit){
  RhatData <- c( 
                 #summary(fit, par="beta_rt")$summary[,"Rhat"],
                 summary(fit, par="beta_rk")$summary[,"Rhat"],
                 summary(fit, par="pai_0r")$summary[,"Rhat"],
                 summary(fit, par="pai_1r")$summary[,"Rhat"],
                 summary(fit, par="theta")$summary[,"Rhat"])
  return(list(meanRhat = mean(RhatData), maxRhat = max(RhatData), countOver11 = sum(RhatData > 1.1)))
}

generate_true_param <- function(setting){
  theta <- rnorm(setting$n_person, 0, 1.0)
  #beta_rt <- generate_beta_rt_self(setting$n_rater, setting$n_time)
  beta_rk <- gen_category_param_previous(setting$K)
  pai_0r <- gen_pai(setting$n_rater)
  #pai_0r <-rnorm(setting$n_rater,0,1)
  pai_1r <- rnorm(setting$n_rater,0,1)
  theta = theta - mean(theta)
  param = list(theta = theta,  
               #beta_rt = beta_rt,
               beta_rk = beta_rk, 
               pai_1r = pai_1r,
               pai_0r = pai_0r
               )
  return(param)
}

get_error <- function(true_param, est_param){
  #error_beta_rt <- list(RMSE = sqrt(apply((true_param$beta_rt-est_param$beta_rt)^2, 2, mean)),
                         #BIAS = apply((true_param$beta_rt - est_param$beta_rt), 2, mean));
  error_pai_0r <- calc_error(true_param$pai_0r, est_param$pai_0r)
  error_pai_1r <- calc_error(true_param$pai_1r, est_param$pai_1r)
  error_category <- list(RMSE = sqrt(apply((true_param$beta_rk - est_param$beta_rk)^2, 2, mean)), 
                         BIAS = apply((true_param$beta_rk - est_param$beta_rk), 2, mean));
  error_theta <- calc_error(true_param$theta, est_param$theta)
  rmse <- list(theta = error_theta, 
               #beta_rt = error_beta_rt, 
               pai_0r = error_pai_0r, 
               pai_1r = error_pai_1r, 
               beta_rk = error_category)  
  return(rmse)
}

