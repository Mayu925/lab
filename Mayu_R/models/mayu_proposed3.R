get_prm_list <- function(param, i, r, t){
  return(list(alpha_r = param$alpha_r[r], pai_1r = param$pai_1r[r], pai_0r = param$pai_0r[r],  category_prm = append(0, param$beta_rk[r,])))
}

logit <- function(setting, param, category_prm, x, t){
  return(param$alpha_r *(x - param$pai_0r - param$pai_1r * (t /setting$n_time)- category_prm))
}

get_param_size <- function(data){
  return(3* data$R + data$R * (data$K - 2) + data$T * 2  +data$I * 2 + data$J - (5 + data$T))
}

get_estimates <- function(fit, setting){
  theta <- summary(fit, par="theta")$summary[,"mean"]
  alpha_r <- get_alpha_estimates_with_restriction(summary(fit, par="alpha_r")$summary[,"mean"])
  category_prm <- convert_category_estimates(summary(fit, par="beta_rk")$summary[,"mean"],setting$n_rater, setting$K)
  pai_0r <- summary(fit, par="pai_0r")$summary[,"mean"]
  pai_1r <- get_alpha_estimates_with_restriction(summary(fit, par="pai_1r")$summary[,"mean"])
  param = list(theta = theta, 
               alpha_r = alpha_r, 
               beta_rk = category_prm, 
               pai_0r = pai_0r, 
               pai_1r = pai_1r)
  return(param)
}

get_Rhat_stat <- function(fit){
  RhatData <- c( 
                 summary(fit, par="alpha_r")$summary[,"Rhat"],
                 summary(fit, par="beta_rk")$summary[,"Rhat"],
                 summary(fit, par="pai_0r")$summary[,"Rhat"],
                 summary(fit, par="pai_1r")$summary[,"Rhat"],
                 summary(fit, par="theta")$summary[,"Rhat"])
  return(list(meanRhat = mean(RhatData), maxRhat = max(RhatData), countOver11 = sum(RhatData > 1.1)))
}

generate_true_param <- function(setting){
  theta <- rnorm(setting$n_person, 0, 1.0)
  alpha_r <- generate_constrained_alpha(setting$n_rater)
  beta_rk <- gen_category_param_previous(setting$K)
  pai_0r <- rnorm(setting$n_rater, 0, 1.0)
  pai_1r <- generate_constrained_alpha(setting$n_rater)
  theta = theta - mean(theta)
  pai_0r = pai_0r - mean(pai_0r)
  param = list(theta = theta,  
               alpha_r = alpha_r, 
               beta_rk = beta_rk, 
               pai_0r = pai_0r, 
               pai_1r = pai_1r)
  return(param)
}

get_error <- function(true_param, est_param){
  error_alpha_r <- calc_error(true_param$alpha_r, est_param$alpha_r)
  error_pai_0r <- calc_error(true_param$pai_0r, est_param$pai_0r)
  error_pai_1r <- calc_error(true_param$pai_1r, est_param$pai_1r)
  error_category <- list(RMSE = sqrt(apply((true_param$beta_rk - est_param$beta_rk)^2, 2, mean)), 
                         BIAS = apply((true_param$beta_rk - est_param$beta_rk), 2, mean));
  error_theta <- calc_error(true_param$theta, est_param$theta)
  rmse <- list(theta = error_theta, 
               alpha_r = error_alpha_r, 
               pai_0r = error_pai_0r, 
               pai_1r = error_pai_1r, 
               beta_rk = error_category)  
  return(rmse)
}

