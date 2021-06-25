get_result_statistics_common <- function(fit1, data1, setting){
  est_param <- get_estimates(fit1, setting)
  Rhat <- get_Rhat_stat(fit1)
  SE_stat <- get_theta_se_stat(fit1)
  D <- c("WAIC", waic(extract_log_lik(fit1))$estimates[3,1], 
         "ML", get_ml(fit1),
         "AIC", aic(data1, est_param),
         "BIC", bic(data1, est_param),
         "RhatMax", Rhat$maxRhat,
         "RhatMean", Rhat$meanRhat,
         "SE Mean", SE_stat$mean_se,
         "SE SD", SE_stat$sd_se
  )
  return(D)
}
