library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("models/ctl_util.R")

model = "mayu_proposed"

source(paste("models/", "mayu_proposed", ".R", sep=""))
stan <- stan_model(file=paste("stan/", model, ".stan", sep=""))
k=5
loop=1
j=100
t=5
r=15
i=1

for(loop in 1:10){
  TH <- c()
  for(j in c(60, 90, 120)){
      for(r in c(10, 15)){
        for (t in c(3,5,10)){
            print(paste(loop, j, 1, r, t, k, sep=","))
            setting <- list(K = k, n_person = j, n_item = 1, n_rater = r, n_time = t)
            true_param <-generate_true_param(setting)
            data <- generate_data(true_param, setting)
            fit <- sampling(stan, data=data, iter=3000, warmup=2000, chains=3, seed=1)
            est_param <- get_estimates(fit, setting)
            d <-  get_error(true_param, est_param)
            Rhat <- get_Rhat_stat(fit)
            TH <- rbind(TH, c(j, r, t,
                            d$theta$RMSE, 
                            d$alpha_r$RMSE, 
                            mean(d$beta_rt$RMSE), 
                            mean(d$beta_rk$RMSE),  
                            d$theta$BIAS, 
                            d$alpha_r$BIAS, 
                            mean(d$beta_rt$BIAS),
                            mean(d$beta_rk$BIAS),
                            Rhat$meanRhat, Rhat$maxRhat))
      }
    }    
  }
  write.csv(TH, paste("output/parameter_recovery/mayu/", model, "/loop_", loop, ".csv", sep=""), row.names = FALSE)
}

TH[1,]
TH <-as.matrix(read.csv(paste("output/parameter_recovery/mayu/", model, "/loop_", 1, ".csv", sep="")))

for(loop in 1:5){

  TH <-TH + as.matrix(read.csv(paste("output/parameter_recovery/mayu/", model, "/loop_", loop, ".csv", sep="")))
}
TH <- TH / 6
write.csv(TH, paste("output/parameter_recovery/mayu/", model, "/parameter_recovery_summary.csv", sep=""), row.names = FALSE)

# 結果のプロット用関数群
plot(true_param$theta, est_param$theta, main="theta")
abline(coef = c(0,1))

plot(true_param$alpha_r, est_param$alpha_r, main="alpha_r")
abline(coef = c(0,1))

plot(true_param$beta_rk[,1], est_param$beta_rk[,1], main="beta_r1")
abline(coef = c(0,1))

plot(true_param$beta_rk[,2], est_param$beta_rk[,2], main="beta_r2")
abline(coef = c(0,1))

plot(true_param$beta_rk[,3], est_param$beta_rk[,3], main="beta_r3")
abline(coef = c(0,1))

plot(true_param$beta_rk[,3], est_param$beta_rk[,4], main="beta_r4")
abline(coef = c(0,1))


plot_beta_rt <- function(r){
  ylim <- c(min(true_param$beta_rt[r,])-0.5, max(true_param$beta_rt[r,])+0.5)
  plot(true_param$beta_rt[r,], xlab="TimeID", ylab="beta_rt" ,type="l",xlim=c(1, 5), xaxp=c(1,5, 4), ylim=ylim,  main=paste("beta_rt(r = ", r, ")", sep=""))
  par(new=T)
  plot(est_param$beta_rt[r,],xlab="TimeID", ylab="beta_rt",xlim=c(1, 5), xaxp=c(1, 5, 4),type="l", ylim=ylim, lty=2,)
}

for(r in 1:15){
  plot_beta_rt(r)
}

def <- list(xlim = 2, xlab = expression(paste("Ability ", theta)),
            ylim = 1.2, ylab = expression(paste("Probability ", P[irk])),
            title = "", SeYlim = 4, legendNcol = 5, legendPosition =
              "topleft", maincx = 1.3 , lcex = 1.3 , axcx = 1.2 , llwd = 3,
            style = c(1, 1, 1, 1, 1), 
            color = c(1, 2, 3, 4, 5),
            pchs = c(0, 1, 2, 3, 4, 5), 
            legends = c("k=1","k=2","k=3","k=4","k=5"))

par(mfrow=c(1,1))
par(mar = c(3, 4.0, 2.1, 0.5))
par(mgp = c(2.2, 0.8, 0))
param1 <- list(alpha_r = 1.0, 
              beta_rt = 0.0, 
              category_prm = c(0,-1.0,0.0,0.5,1.0))
param2 <- list(alpha_r = 1.0, 
              beta_rt = 0.5, 
              category_prm = c(0,-1.0,0.0,0.5,1.0))
param3 <- list(alpha_r = 1.0, 
              beta_rt = -0.5, 
              category_prm = c(0,-1.0,0.0,0.5,1.0))
param4 <- list(alpha_r = 1.0, 
              beta_rt = 0.0, 
              category_prm = c(0,-1.0,0.0,0.2,1.0))
param5 <- list(alpha_r = 2.0, 
              beta_rt = 0.0, 
              category_prm = c(0,-1.0,0.0,0.5,1.0))
param6 <- list(alpha_r = 0.5, 
              beta_rt = 0.0, 
              category_prm = c(0,-1.0,0.0,0.5,1.0))
draw_icc(param1, 5, def)
legend(def$legendPosition, cex = def$axcx, lwd=def$llwd,
       legend=def$legends, col=def$color, lty=def$style, bg = "white", ncol =
         def$legendNcol, pch=def$pchs)
draw_icc(param2, 5, def)
legend(def$legendPosition, cex = def$axcx, lwd=def$llwd,
       legend=def$legends, col=def$color, lty=def$style, bg = "white", ncol =
         def$legendNcol, pch=def$pchs)
draw_icc(param3, 5, def)
legend(def$legendPosition, cex = def$axcx, lwd=def$llwd,
       legend=def$legends, col=def$color, lty=def$style, bg = "white", ncol =
         def$legendNcol, pch=def$pchs)
draw_icc(param4, 5, def)
legend(def$legendPosition, cex = def$axcx, lwd=def$llwd,
       legend=def$legends, col=def$color, lty=def$style, bg = "white", ncol =
         def$legendNcol, pch=def$pchs)
draw_icc(param5, 5, def)
legend(def$legendPosition, cex = def$axcx, lwd=def$llwd,
       legend=def$legends, col=def$color, lty=def$style, bg = "white", ncol =
         def$legendNcol, pch=def$pchs)
draw_icc(param6, 5, def)
legend(def$legendPosition, cex = def$axcx, lwd=def$llwd,
       legend=def$legends, col=def$color, lty=def$style, bg = "white", ncol =
         def$legendNcol, pch=def$pchs)