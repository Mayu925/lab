#Sys.setenv("http_proxy" = "http://proxy.uec.ac.jp:8080")
#Sys.setenv("https_proxy" = "https://proxy.uec.ac.jp:8080")
library(rstan)
library(loo)
library(psych)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("models/ctl_util_previous.R")

model = "mayu_previous"

setting1 <- list(K = 5, n_person = 134, n_rater = 17, n_time = 4)
setting2 <- list(K = 5, n_person = 34, n_item = 1, n_rater = 34, n_time = 5)
setting3 <- list(K = 5, n_person = 34, n_item = 1, n_rater = 34, n_time = 10)
data1 <- read_data(setting1, paste("data/mayu_data_3_t=3.csv", sep=""))
data2 <-read_data(setting2, paste("data/mayu_data_3_t=5.csv", sep=""))
data3 <-read_data(setting3, paste("data/mayu_data_3_t=10.csv", sep=""))
data1 <-read_data(setting1, paste("data/data_o_2-5-11.csv", sep=""))
data1$ItemID = data1$ItemID -2
data2$ItemID = data2$ItemID -2
data3$ItemID = data3$ItemID -2

stan <- stan_model(file=paste("stan/", model, ".stan", sep=""))
fit1 <- sampling(stan, data=data1, iter=1000, warmup=500, chains=3)
fit2 <- sampling(stan, data=data2, iter=1000, warmup=500, chains=3)
fit3 <- sampling(stan, data=data3, iter=1000, warmup=500, chains=3)
#fit4 <- sampling(stan, data=data4, iter=1000, warmup=500, chains=3)

source(paste( "models/", "mayu_previous", ".R", sep=""))
est_param1 <- get_estimates(fit1, setting1)
D1 <- get_result_statistics_common(fit1, data1, setting1)
write.csv(t(matrix(D1, nrow=2)), paste("output/realdata/MCMC_mayu/previous/", model, "t_3.csv", sep=""))
#write.csv(t(matrix(est_param1, nrow=34)), paste("output/realdata/parameters/mayu/previous/", model, "t_3.csv", sep=""))

est_param2 <- get_estimates(fit2, setting2)
D2 <- get_result_statistics_common(fit2, data2, setting2)
write.csv(t(matrix(D2, nrow=2)), paste("output/realdata/MCMC_mayu/previous/", model, "t_5.csv", sep=""))
#write.csv(t(matrix(est_param2$theta, nrow=34)), paste(dir, "output/realdata/parameters/mayu/previous/", model, "t_5.csv", sep=""))

est_param3 <- get_estimates(fit3, setting3)
D3 <- get_result_statistics_common(fit3, data3, setting3)
write.csv(t(matrix(D3, nrow=2)), paste("output/realdata/MCMC_mayu/previous/", model, "t_10.csv", sep=""))
#write.csv(t(matrix(est_param1$theta, nrow=34)), paste(dir, "output/realdata/parameters/mayu/previous", model, "t_10.csv", sep=""))

wbic3 <- -mean(rowSums(extract(fit1)$log_lik))
wbic5 <- -mean(rowSums(extract(fit2)$log_lik))
wbic10 <- -mean(rowSums(extract(fit3)$log_lik))
wbic3
wbic5
wbic10
SD <- summary(fit1, par="theta")$summary[,c("sd", "se_mean")]
apply(SD, 2, mean)
write.csv(SD, paste(dir, "output/realdata/SE/", model, "1.csv", sep=""))

sink(paste(dir, "output/realdata/parameters/", model, "_param.txt", sep=""))
sink()

