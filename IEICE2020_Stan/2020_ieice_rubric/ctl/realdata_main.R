library(rstan)
library(loo)
library(psych)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dir <- "2020_ieice_rubric/"

model = "gpcm_rubric_uto"

source("common/util.R")
source("common/ctl_util.R")

setting <- list(K = 4, n_person = 34, n_item = 4, n_rater = 10, n_rubric = 5)
data <- read_data(setting, paste("data/rubric_data/rubric.csv", sep=""))

stan <- stan_model(file=paste(dir, "stan/", model, ".stan", sep=""))
fit <- sampling(stan, data=data, iter=1000, warmup=500, chains=3)

source(paste(dir, "models/", model, ".R", sep=""))
est_param <- get_estimates(fit, setting)
D <- get_result_statistics_common(fit, data, setting)
write.csv(t(matrix(D, nrow=2)), paste(dir, "output/realdata/MCMC_statistics/", model, ".csv", sep=""))

SD <- summary(fit, par="theta")$summary[,c("sd", "se_mean")]
apply(SD, 2, mean)
write.csv(SD, paste(dir, "output/realdata/SE/", model, ".csv", sep=""))

sink(paste(dir, "output/realdata/parameters/", model, "_param.txt", sep=""))
est_param
sink()
