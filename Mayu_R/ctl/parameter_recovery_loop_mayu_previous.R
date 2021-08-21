library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dir <- ""

source("models/ctl_util.R")

model = "mayu_previous"

source(paste(dir, "models/", model, ".R", sep=""))
stan <- stan_model(file=paste(dir, "stan/", model, ".stan", sep=""))
k=5
loop=1
j=20
r=5

for(loop in 1:2){
  TH <- c()
  for(j in c(30, 50)){
    for(r in c(5, 10)){
            print(paste(loop, j, 1, r, j, k, sep=","))
            setting <- list(K = k, n_person = j, n_item = 1, n_rater = r, n_time = j)
            true_param <-generate_true_param(setting)
            data <- generate_data(true_param, setting)
            fit <- sampling(stan, data=data, iter=1000, warmup=500, chains=3, seed=1)
            est_param <- get_estimates(fit, setting)
            d <-  get_error(true_param, est_param)
            Rhat <- get_Rhat_stat(fit)
            TH <- rbind(TH, c(j, r, k, 
                            d$theta$RMSE, 
                            d$pai_0r$RMSE,
                            d$pai_1r$RMSE,
                            #mean(d$alpha_rt$RMSE), 
                            mean(d$category_prm$RMSE),  
                            d$theta$BIAS, 
                            d$pai_0r$BIAS,
                            d$pai_1r$BIAS,
                            #mean(d$alpha_rt$BIAS),
                            mean(d$category_prm$BIAS),
                            Rhat$meanRhat, Rhat$maxRhat))
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


