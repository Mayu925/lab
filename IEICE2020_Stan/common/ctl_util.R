get_result_statistics_common <- function(fit, data, setting){
  est_param <- get_estimates(fit, setting)
  Rhat <- get_Rhat_stat(fit)
  SE_stat <- get_theta_se_stat(fit)
  D <- c("WAIC", waic(extract_log_lik(fit))$estimates[3,1], 
         "ML", get_ml(fit),
         "AIC", aic(data, est_param),
         "BIC", bic(data, est_param),
         "RhatMax", Rhat$maxRhat,
         "RhatMean", Rhat$meanRhat,
         "SE Mean", SE_stat$mean_se,
         "SE SD", SE_stat$sd_se
  )
  return(D)
}
