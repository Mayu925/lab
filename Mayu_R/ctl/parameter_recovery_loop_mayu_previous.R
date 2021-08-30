library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("models/ctl_util.R")

model = "mayu_previous"

source(paste("models/", model, ".R", sep=""))
stan <- stan_model(file=paste("stan/", model, ".stan", sep=""))
loop=1
j=20
i=1
r=5
t=j
k=5
for(loop in 1:3){
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
                            mean(d$alpha_rt$RMSE), 
                            mean(d$beta_rk$RMSE),  
                            d$theta$BIAS, 
                            d$pai_0r$BIAS,
                            d$pai_1r$BIAS,
                            mean(d$alpha_rt$BIAS),
                            mean(d$beta_rk$BIAS),
                            Rhat$meanRhat, Rhat$maxRhat))
    }
  }    
  write.csv(TH, paste("output/parameter_recovery/mayu/", model, "/loop_", loop, ".csv", sep=""), row.names = FALSE)
}

TH[1,]
TH <-as.matrix(read.csv(paste("output/parameter_recovery/mayu/", model, "/loop_", loop, ".csv", sep="")))
for(loop in 1:4){
  TH <-TH + as.matrix(read.csv(paste("output/parameter_recovery/mayu/", model, "/loop_", loop, ".csv", sep="")))
}
TH <- TH / 5
write.csv(TH, paste("output/parameter_recovery/mayu/", model, "/parameter_recovery_summary.csv", sep=""), row.names = FALSE)

# 結果のプロット用関数群
plot(true_param$theta, est_param$theta, main="theta")
abline(coef = c(0,1))

plot(true_param$pai_0r, est_param$pai_0r, main="pai_0r")
abline(coef = c(0,1))

plot(true_param$pai_1r, est_param$pai_1r, main="pai_1r")
abline(coef = c(0,1))

plot_alpha_rt <- function(r){
  ylim <- c(min(true_param$alpha_rt[r,])-0.5, max(true_param$alpha_rt[r,])+0.5)
  plot(true_param$alpha_rt[r,], type="l", ylim=ylim, ylab="", main=paste("alpha_rt(r = ", r, ")", sep=""))
  par(new=T)
  plot(est_param$alpha_rt[r,], type="l", ylim=ylim, lty=2, ylab="")
}
for(r in 1:10){
  plot_alpha_rt(r)
}

plot(true_param$beta_rk, est_param$beta_rk, main="beta_rk")
abline(coef = c(0,1))
