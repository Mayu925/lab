library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dir <- "2020_ieice_rubric/"

source("common/util_mayu.R")

model = "mayu_new"

source(paste(dir, "models/", model, ".R", sep=""))
stan <- stan_model(file=paste(dir, "stan/", model, ".stan", sep=""))

for(loop in 28:30){
  TH <- c()
  for(j in c(30, 50)){
    for(i in c(3, 5)){
      for(r in c(5, 10)){
        for(t in c(5, 10)){
          for(k in c(3, 4)){
            print(paste(loop, j, i, r, t, k, sep=","))
            setting <- list(K = k, n_person = j, n_item = i, n_rater = r, n_time = t)
            true_param <-generate_true_param(setting)
            data <- generate_data(true_param, setting)
            fit <- sampling(stan, data=data, iter=1000, warmup=500, chains=3, seed=1)
            est_param <- get_estimates(fit, setting)
            d <-  get_error(true_param, est_param)
            Rhat <- get_Rhat_stat(fit)
            TH <- rbind(TH, c(j, i, r, t, k, 
                            d$theta$RMSE, 
                            d$alpha_r$RMSE, 
                            #mean(d$alpha_rt$RMSE), 
                            mean(d$beta_rk$RMSE),  
                            d$theta$BIAS, 
                            d$alpha_r$BIAS, 
                            #mean(d$alpha_rt$BIAS),
                            mean(d$beta_rk$BIAS),
                            Rhat$meanRhat, Rhat$maxRhat))
          }
        }
      }
    }
  }    
  write.csv(TH, paste(dir, "output/parameter_recovery/mayu/", model, "/loop_", loop, ".csv", sep=""), row.names = FALSE)
}

TH[1,]
TH <-as.matrix(read.csv(paste(dir, "output/parameter_recovery/mayu/", model, "/loop_", 28, ".csv", sep="")))
for(loop in 29:30){
  TH <-TH + as.matrix(read.csv(paste(dir, "output/parameter_recovery/mayu/", model, "/loop_", loop, ".csv", sep="")))
}
TH <- TH / 3
write.csv(TH, paste(dir, "output/parameter_recovery/mayu/", model, "/parameter_recovery_summary.csv", sep=""), row.names = FALSE)


