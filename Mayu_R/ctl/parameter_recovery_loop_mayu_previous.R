library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("models/ctl_util_previous.R")

model = "mayu_previous"

source(paste("models/", model, ".R", sep=""))

stan <- stan_model(file=paste("stan/", model, ".stan", sep=""))
loop=1
j=30
i=1
r=10
t=j
k=5

for(loop in 1:10){
  TH <- c()
  for(j in c(60, 90, 120)){
    for (t in c(3,5,10)){
    for(r in c(10, 15)){
            print(paste(loop, j, 1, r, t, k, sep=","))
            setting <- list(K = k, n_person = j, n_rater = r, n_time = t)
            true_param <-generate_true_param(setting)
            data <- generate_data(true_param, setting)
            fit <- sampling(stan, data=data, iter=1000, warmup=500, chains=3, seed=1)
            est_param <- get_estimates(fit, setting)
            d <-  get_error(true_param, est_param)
            Rhat <- get_Rhat_stat(fit)
            TH <- rbind(TH, c(j, r, t, 
                            d$theta$RMSE, 
                            d$pai_0r$RMSE,
                            d$pai_1r$RMSE,
                            #mean(d$beta_rt$RMSE), 
                            mean(d$beta_rk$RMSE),  
                            d$theta$BIAS, 
                            d$pai_0r$BIAS,
                            d$pai_1r$BIAS,
                            #mean(d$beta_rt$BIAS),
                            mean(d$beta_rk$BIAS),
                            Rhat$meanRhat, Rhat$maxRhat))
    }
    }    
  }
  write.csv(TH, paste("output/parameter_recovery/mayu/", model, "/loop2_", loop, ".csv", sep=""), row.names = FALSE)
}

TH[1,]
TH <-as.matrix(read.csv(paste("output/parameter_recovery/mayu/", model, "/loop2_", 1, ".csv", sep="")))
for(loop in 2:10){
  TH <-TH + as.matrix(read.csv(paste("output/parameter_recovery/mayu/", model, "/loop2_", loop, ".csv", sep="")))
}
#TH <-TH + as.matrix(read.csv(paste("output/parameter_recovery/mayu/", model, "/loop_", 1, ".csv", sep="")))
TH <- TH / 10
write.csv(TH, paste("output/parameter_recovery/mayu/", model, "/parameter_recovery_summary2.csv", sep=""), row.names = FALSE)

# 結果のプロット用関数群
plot(true_param$theta, est_param$theta, main="theta")
abline(coef = c(0,1))

plot(true_param$pai_0r, est_param$pai_0r, main="pai_0r")
abline(coef = c(0,1))

plot(true_param$pai_1r, est_param$pai_1r, main="pai_1r")
abline(coef = c(0,1))

plot_beta_rt <- function(r){
  ylim <- c(min(true_param$beta_rt[r,])-0.5, max(true_param$beta_rt[r,])+0.5)
  plot(true_param$beta_rt[r,], type="l", ylim=ylim, ylab="", main=paste("beta_rt(r = ", r, ")", sep=""))
  par(new=T)
  plot(est_param$beta_rt[r,], type="l", ylim=ylim, lty=2, ylab="")
}
for(r in 1:10){
  plot_beta_rt(r)
}
plot(true_param$beta_rk, est_param$beta_rk, main="beta_rk")
abline(coef = c(0,1))
