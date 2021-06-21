#Sys.setenv("http_proxy" = "http://proxy.uec.ac.jp:8080")
#Sys.setenv("https_proxy" = "https://proxy.uec.ac.jp:8080")
library(rstan)
library(loo)
library(psych)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dir <- "2020_ieice_rubric/"

model = "mayu"

source("common/util_mayu.R")
source("common/ctl_util.R")


setting <- list(K = 5, n_person = 34, n_item = 4, n_rater = 34, n_time = 34)
data <- read_data(setting, paste("data/mayu_data_0.csv", sep=""))
data1 <-read_data(setting, paste("data/mayu_data_1.csv", sep=""))
data2 <-read_data(setting, paste("data/mayu_data_2.csv", sep=""))
data3 <-read_data(setting, paste("data/mayu_data_3.csv", sep=""))
data$ItemID = data$ItemID + 1
data1$TimeID = data1$TimeID - 34
data2$TimeID = data1$TimeID - 68
data3$TimeID = data1$TimeID - 102
data1<- list()
data2<- list()
data3<- list()
data4<- list()
data1$K=data$K
data2$K=data$K
data3$K=data$K
data4$K=data$K
data1$J=data$J
data2$J=data$J
data3$J=data$J
data4$J=data$J
data1$I=4
data2$I=4
data3$I=4
data4$I=4
data1$R=data$R
data2$R=data$R
data3$R=data$R
data4$R=data$R
data1$T=data$T
data2$T=data$T
data3$T=data$T
data4$T=data$T
data1$N=1156
data2$N=1156
data3$N=1156
data4$N=1156

j1=1
j2=1
j3=1
j4=1

for(i in 1:4624){
  if(data$ItemID[i] == 1){
    data1$ItemID[j1]=data$ItemID[i]
    data1$ExamineeID[j1]=data$ExamineeID[i]
    data1$RaterID[j1]=data$RaterID[i]
    data1$TimeID[j1]=data$TimeID[i]
    data1$X[j1]=data$X[i]
    j1=j1+1
  }else if(data$ItemID[i] == 2){
    data2$ItemID[j2]=data$ItemID[i]
    data2$ExamineeID[j2]=data$ExamineeID[i]
    data2$RaterID[j2]=data$RaterID[i]
    data2$TimeID[j2]=data$TimeID[i]
    data2$X[j2]=data$X[i]
    j2=j2+1
  }else if(data$ItemID[i] == 3){
    data3$ItemID[j3]=data$ItemID[i]
    data3$ExamineeID[j3]=data$ExamineeID[i]
    data3$RaterID[j3]=data$RaterID[i]
    data3$TimeID[j3]=data$TimeID[i]
    data3$X[j3]=data$X[i]
    j3=j3+1
  }else if(data$ItemID[i] == 4){
    data4$ItemID[j4]=data$ItemID[i]
    data4$ExamineeID[j4]=data$ExamineeID[i]
    data4$RaterID[j4]=data$RaterID[i]
    data4$TimeID[j4]=data$TimeID[i]
    data4$X[j4]=data$X[i]
    j4=j4+1
  }
}


stan <- stan_model(file=paste(dir, "stan/", model, ".stan", sep=""))
fit1 <- sampling(stan, data=data, iter=1000, warmup=500, chains=3)
summary(fit1)$summary[,"mean"]
fit2 <- sampling(stan, data=data2, iter=1000, warmup=500, chains=3)
fit3 <- sampling(stan, data=data3, iter=1000, warmup=500, chains=3)
fit4 <- sampling(stan, data=data4, iter=1000, warmup=500, chains=3)

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
