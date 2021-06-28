Sys.setenv("http_proxy" = "http://proxy.uec.ac.jp:8080")
Sys.setenv("https_proxy" = "https://proxy.uec.ac.jp:8080")
library(rstan)
library(loo)
library(psych)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dir <- "2020_ieice_rubric/"

#model = "mayu_new"
model = "mayu"

source("common/util_mayu.R")
source("common/ctl_util.R")


setting <- list(K = 5, n_person = 34, n_item = 4, n_rater = 34, n_time = 34)
data1 <- read_data(setting, paste("data/mayu_data_0.csv", sep=""))
data2 <-read_data(setting, paste("data/mayu_data_1.csv", sep=""))
data3 <-read_data(setting, paste("data/mayu_data_2.csv", sep=""))
data4 <-read_data(setting, paste("data/mayu_data_3.csv", sep=""))

data1$ItemID = data1$ItemID + 1
data2$ItemID = data2$ItemID + 1
data3$ItemID = data3$ItemID + 1
data4$ItemID = data4$ItemID + 1
data2$TimeID = data2$TimeID - 34
data3$TimeID = data3$TimeID - 68
data4$TimeID = data4$TimeID - 102


stan <- stan_model(file=paste(dir, "stan/", model, ".stan", sep=""))
fit1 <- sampling(stan, data=data1, iter=1000, warmup=500, chains=3)
summary(fit1)$summary[,"mean"]
fit2 <- sampling(stan, data=data2, iter=1000, warmup=500, chains=3)
fit3 <- sampling(stan, data=data3, iter=1000, warmup=500, chains=3)
fit4 <- sampling(stan, data=data4, iter=1000, warmup=500, chains=3)

source(paste(dir, "models/", model, ".R", sep=""))
est_param1 <- get_estimates(fit1, setting)
D1 <- get_result_statistics_common(fit1, data1, setting)
write.csv(t(matrix(D1, nrow=2)), paste(dir, "output/realdata/MCMC_statistics/", model, "1.csv", sep=""))

est_param2 <- get_estimates(fit2, setting)
D2 <- get_result_statistics_common(fit2, data2, setting)
write.csv(t(matrix(D2, nrow=2)), paste(dir, "output/realdata/MCMC_statistics/", model, "2.csv", sep=""))

est_param3 <- get_estimates(fit3, setting)
D3 <- get_result_statistics_common(fit3, data3, setting)
write.csv(t(matrix(D3, nrow=2)), paste(dir, "output/realdata/MCMC_statistics/", model, "3.csv", sep=""))

est_param4 <- get_estimates(fit4, setting)
D4 <- get_result_statistics_common(fit4, data4, setting)
write.csv(t(matrix(D4, nrow=2)), paste(dir, "output/realdata/MCMC_statistics/", model, "4.csv", sep=""))

SD <- summary(fit, par="theta")$summary[,c("sd", "se_mean")]
apply(SD, 2, mean)
write.csv(SD, paste(dir, "output/realdata/SE/", model, ".csv", sep=""))

sink(paste(dir, "output/realdata/parameters/", model, "_param.txt", sep=""))
est_param
sink()
