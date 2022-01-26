#Sys.setenv("http_proxy" = "http://proxy.uec.ac.jp:8080")
#Sys.setenv("https_proxy" = "https://proxy.uec.ac.jp:8080")
library(rstan)
library(loo)
library(psych)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

model = "mayu_proposed"

source("models/ctl_util.R")

setting <- list(K = 5, n_person = 134, n_rater = 6, n_time = 4)
data <-read_data(setting, paste("data/data_bias.csv", sep=""))

setting <- list(K = 5, n_person = 134, n_rater = 11, n_time = 4)
data <-read_data(setting, paste("data/data_n_bias.csv", sep=""))

setting <- list(K = 5, n_person = 134, n_rater = 16, n_time = 4)
data <- read_data(setting, paste("data/data_all.csv", sep=""))


stan <- stan_model(file=paste("stan/", model, ".stan", sep=""))
fit <- sampling(stan, data=data, iter=1000, warmup=500, chains=3)

source(paste("models/", "mayu_proposed", ".R", sep=""))

est_param <- get_estimates(fit, setting)
D <- get_result_statistics_common(fit, data, setting)
write.csv(t(matrix(D, nrow=2)), paste( "output/realdata/MCMC_mayu/proposed/", model, ".csv", sep=""))

write.csv(t(matrix(est_param$theta, nrow=134)), paste( "output/realdata/parameters/mayu/proposed/", model, "_theta.csv", sep=""))
write.csv(t(matrix(est_param$alpha_r, nrow=16)), paste( "output/realdata/parameters/mayu/proposed/", model, "_alpha_r.csv", sep=""))
write.csv(t(matrix(est_param$beta_rt, nrow=16)), paste( "output/realdata/parameters/mayu/proposed/", model, "_beta_rt.csv", sep=""))
write.csv(t(matrix(est_param$beta_rk, nrow=16)), paste( "output/realdata/parameters/mayu/proposed/", model, "_beta_rk.csv", sep=""))

wbic <- -mean(rowSums(extract(fit)$log_lik))
wbic

SD <- summary(fit, par="theta")$summary[,c("sd", "se_mean")]
apply(SD, 2, mean)
write.csv(SD, paste(dir, "output/realdata/SE_mayu/", model, "1.csv", sep=""))

sink(paste(dir, "output/realdata/parameters/mayu/proposed/", model, "_param.txt", sep=""))
sink()

# 結果描画用
plot(est_param$theta)
plot(est_param$alpha_r)
est_param$beta_rt
for(r in 1:10){
  plot(est_param$beta_rt[r,],xlab="TimeID",ylab="beta_rt", cex.axis = 1.5,cex.lab=1.8,type="l",xlim=c(1, 4), xaxp=c(1, 4, 3), ylim=c(-1.2,0.8), col=rgb(0.7,0.7,0.7), lwd = 2)
  par(new=T)
  }

for(r in 11){
  plot(est_param$beta_rt[r,],xlab="TimeID",ylab="beta_rt", cex.axis = 1.5,cex.lab=1.8,type="l",xlim=c(1, 4), xaxp=c(1, 4, 3), ylim=c(-1.2,0.8), col=2, lwd = 2)
  par(new=T)
}
for(r in 12:16){
  plot(est_param$beta_rt[r,],xlab="TimeID",ylab="beta_rt", cex.axis = 1.5,cex.lab=1.8,type="l",xlim=c(1, 4), xaxp=c(1, 4, 3), ylim=c(-1.2,0.8), col=rgb(0.7,0.7,0.7), lwd = 2)
  par(new=T)
}
