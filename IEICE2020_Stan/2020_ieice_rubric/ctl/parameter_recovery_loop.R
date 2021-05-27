library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dir <- "2020_ieice_rubric/"

source("common/util.R")

model = "gpcm_rubric_uto"

source(paste(dir, "models/", model, ".R", sep=""))
stan <- stan_model(file=paste(dir, "stan/", model, ".stan", sep=""))

for(loop in 28:30){
  TH <- c()
  for(j in c(30, 50)){
    for(i in c(3, 5)){
      for(r in c(5, 10)){
        for(c in c(5, 10)){
          for(k in c(3, 4)){
            print(paste(loop, j, i, r, c, k, sep=","))
            setting <- list(K = k, n_person = j, n_item = i, n_rater = r, n_rubric = c)
            true_param <-generate_true_param(setting)
            data <- generate_data(true_param, setting)
            fit <- sampling(stan, data=data, iter=1000, warmup=500, chains=3, seed=1)
            est_param <- get_estimates(fit, setting)
            d <-  get_error(true_param, est_param)
            Rhat <- get_Rhat_stat(fit)
            TH <- rbind(TH, c(j, i, r, c, k, 
                            d$theta$RMSE, 
                            d$alpha_i$RMSE, 
                            d$alpha_r$RMSE, 
                            d$alpha_c$RMSE, 
                            d$beta_i$RMSE, 
                            d$beta_r$RMSE, 
                            d$beta_c$RMSE, 
                            d$tau_r$RMSE, 
                            mean(d$beta_ck$RMSE), 
                            d$theta$BIAS, 
                            d$alpha_i$BIAS, 
                            d$alpha_r$BIAS, 
                            d$alpha_c$BIAS, 
                            d$beta_i$BIAS, 
                            d$beta_r$BIAS, 
                            d$beta_c$BIAS, 
                            d$tau_r$BIAS, 
                            mean(d$beta_ck$BIAS),
                            Rhat$meanRhat, Rhat$maxRhat))
          }
        }
      }
    }
  }    
  write.csv(TH, paste(dir, "output/parameter_recovery/", model, "/loop_", loop, ".csv", sep=""), row.names = FALSE)
}

TH[1,]
TH <-as.matrix(read.csv(paste(dir, "output/parameter_recovery/", model, "/loop_", 1, ".csv", sep="")))
for(loop in 2:30){
  TH <-TH + as.matrix(read.csv(paste(dir, "output/parameter_recovery/", model, "/loop_", 1, ".csv", sep="")))
}
TH <- TH / 30
write.csv(TH, paste(dir, "output/parameter_recovery/", model, "/parameter_recovery_summary.csv", sep=""), row.names = FALSE)


