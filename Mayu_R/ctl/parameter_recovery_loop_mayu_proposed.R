library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("models/ctl_util.R")

model = "mayu_proposed"

source(paste("models/", model, ".R", sep=""))
stan <- stan_model(file=paste("stan/", model, ".stan", sep=""))
k=5
loop=1
j=30
<<<<<<< HEAD
r=10
=======
r=15
>>>>>>> origin/master

for(loop in 1:2){
  TH <- c()
  for(j in c(30, 50)){
<<<<<<< HEAD
      for(r in c(5, 10)){
=======
      for(r in c(15, 10)){
>>>>>>> origin/master
            print(paste(loop, j, 1, r, j, k, sep=","))
            setting <- list(K = k, n_person = j, n_item = 1, n_rater = r, n_time = j)
            true_param <-generate_true_param(setting)
            data <- generate_data(true_param, setting)
            fit <- sampling(stan, data=data, iter=3000, warmup=2000, chains=3, seed=1)
            est_param <- get_estimates(fit, setting)
            d <-  get_error(true_param, est_param)
            Rhat <- get_Rhat_stat(fit)
            TH <- rbind(TH, c(j, r, k,
                            d$theta$RMSE, 
                            d$alpha_r$RMSE, 
<<<<<<< HEAD
                            mean(d$alpha_rt$RMSE), 
                            mean(d$beta_rk$RMSE),  
                            d$theta$BIAS, 
                            d$alpha_r$BIAS, 
                            mean(d$alpha_rt$BIAS),
=======
                            mean(d$beta_rt$RMSE), 
                            mean(d$beta_rk$RMSE),  
                            d$theta$BIAS, 
                            d$alpha_r$BIAS, 
                            mean(d$beta_rt$BIAS),
>>>>>>> origin/master
                            mean(d$beta_rk$BIAS),
                            Rhat$meanRhat, Rhat$maxRhat))
      }
  }    
  write.csv(TH, paste("output/parameter_recovery/mayu/", model, "/loop_", 1, ".csv", sep=""), row.names = FALSE)
}

TH[1,]
TH <-as.matrix(read.csv(paste("output/parameter_recovery/mayu/", model, "/loop_", 1, ".csv", sep="")))
<<<<<<< HEAD
for(loop in 29:30){
=======
for(loop in 1:2){
>>>>>>> origin/master
  TH <-TH + as.matrix(read.csv(paste("output/parameter_recovery/mayu/", model, "/loop_", loop, ".csv", sep="")))
}
TH <- TH / 3
write.csv(TH, paste("output/parameter_recovery/mayu/", model, "/parameter_recovery_summary.csv", sep=""), row.names = FALSE)

# 結果のプロット用関数群
plot(true_param$theta, est_param$theta, main="theta")
abline(coef = c(0,1))

plot(true_param$alpha_r, est_param$alpha_r, main="alpha_r")
abline(coef = c(0,1))

plot(true_param$category_prm[,1], est_param$category_prm[,1], main="beta_r1")
abline(coef = c(0,1))

plot(true_param$category_prm[,2], est_param$category_prm[,2], main="beta_r2")
abline(coef = c(0,1))

plot(true_param$category_prm[,3], est_param$category_prm[,3], main="beta_r3")
abline(coef = c(0,1))

plot(true_param$category_prm[,3], est_param$category_prm[,4], main="beta_r4")
abline(coef = c(0,1))

<<<<<<< HEAD
plot_alpha_rt <- function(r){
  ylim <- c(min(true_param$alpha_rt[r,])-0.5, max(true_param$alpha_rt[r,])+0.5)
  plot(true_param$alpha_rt[r,], type="l", ylim=ylim, ylab="", main=paste("alpha_rt(r = ", r, ")", sep=""))
  par(new=T)
  plot(est_param$alpha_rt[r,], type="l", ylim=ylim, lty=2, ylab="")
}
for(r in 1:10){
  plot_alpha_rt(r)
=======
plot_beta_rt <- function(r){
  ylim <- c(min(true_param$beta_rt[r,])-0.5, max(true_param$beta_rt[r,])+0.5)
  plot(true_param$beta_rt[r,], type="l", ylim=ylim, ylab="", main=paste("beta_rt(r = ", r, ")", sep=""))
  par(new=T)
  plot(est_param$beta_rt[r,], type="l", ylim=ylim, lty=2, ylab="")
}
for(r in 1:10){
  plot_beta_rt(r)
>>>>>>> origin/master
}
