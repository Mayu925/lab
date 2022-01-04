#Sys.setenv("http_proxy" = "http://proxy.uec.ac.jp:8080")
#Sys.setenv("https_proxy" = "https://proxy.uec.ac.jp:8080")
library(rstan)
library(loo)
library(psych)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("models/ctl_util_previous.R")

model = "mayu_previous_WBIC"


setting <- list(K = 5, n_person = 134, n_rater = 6, n_time = 4)
data <-read_data(setting, paste("data/data_bias.csv", sep=""))

setting <- list(K = 5, n_person = 134, n_rater = 10, n_time = 4)
data <-read_data(setting, paste("data/data_n_bias.csv", sep=""))

setting <- list(K = 5, n_person = 134, n_rater = 16, n_time = 4)
data <- read_data(setting, paste("data/data_all.csv", sep=""))


stan <- stan_model(file=paste("stan/", model, ".stan", sep=""))
fit <- sampling(stan, data=data, iter=1000, warmup=500, chains=3)

source(paste( "models/", "mayu_previous", ".R", sep=""))
est_param <- get_estimates(fit, setting)
D <- get_result_statistics_common(fit, data, setting)
write.csv(t(matrix(D, nrow=2)), paste("output/realdata/MCMC_mayu/previous/", model, ".csv", sep=""))
#write.csv(t(matrix(est_param1, nrow=34)), paste("output/realdata/parameters/mayu/previous/", model, "t_3.csv", sep=""))

wbic <- -mean(rowSums(extract(fit)$log_lik))
wbic

SD <- summary(fit1, par="theta")$summary[,c("sd", "se_mean")]
apply(SD, 2, mean)
write.csv(SD, paste(dir, "output/realdata/SE/", model, "1.csv", sep=""))

sink(paste(dir, "output/realdata/parameters/", model, "_param.txt", sep=""))
sink()
pai_1rt <- matrix(nrow = setting$n_rater, ncol = 4)
for(t in 1:4){
  for(r in 1:setting$n_rater){
    pai_1rt[r,t] = est_param$pai_1r[r]*(t/4)
  }
}
# 結果描画用
plot(est_param$theta)
plot(est_param$beta_rk)
for(r in 1:setting1$n_rater){
  plot(est_param1$pai_0r[r]+pai_1rt[r,],xlab="TimeID",ylab="pai_0r*pai_1r*t", type="l",xlim=c(1, 4), xaxp=c(1, 4, 3), ylim=c(-1,1), col=r, lwd = 1)
  par(new=T)
}
