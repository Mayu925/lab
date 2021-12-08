get_result_statistics_common <- function(fit, data, setting){
  est_param <- get_estimates(fit, setting)
  Rhat <- get_Rhat_stat(fit)
  SE_stat <- get_theta_se_stat(fit)
  D <- c("WAIC", waic(extract_log_lik(fit))$estimates[3,1], 
         "RhatMax", Rhat$maxRhat,
         "RhatMean", Rhat$meanRhat,
         "SE Mean", SE_stat$mean_se,
         "SE SD", SE_stat$sd_se
  )
  return(D)
}

convert_category_estimates <- function(category_prm, N, K){
  for(n in 1:N){
    prm = category_prm[((n-1)*(K-2)+1):((n-1)*(K-2)+(K-2))]
    prm = append(prm, -1*sum(prm))
    if(n == 1){
      mat = t(data.frame(prm))
    } else {
      mat = rbind(mat, t(data.frame(prm)))
    }
  }  
  return(mat)
}

convert_category_estimates_previous <- function(category_prm, K){
    prm = category_prm[1:(K-2)]
    prm = append(prm, -1*sum(prm))
    mat = t(data.frame(prm))
  return(mat)
}

convert_beta_rt <- function(beta_rt, T, R){
  for(r in 1:R){
    alrt = beta_rt[((r-1)*T+1):((r-1)*T+T)]
    if(r == 1){
      mat = t(data.frame(alrt))
    }else{
      mat = rbind(mat,t(data.frame(alrt)))
    }
  }
  return(mat)
}


convert_beta_rt2 <- function(beta_rt,T){
  for(t in 1:T){
    beta_rt[t] = T
  }
  return(beta_rt)
}

prob <-function(param, k, x){
  all_sum <- 0
  tmp <- 0
  for(m in 1:length(param$category_prm)){
    tmp <- tmp + 1.7 * logit(param, param$category_prm[m], x)
    if(m==k) target <- exp(tmp)
    all_sum <- all_sum + exp(tmp)
  }
  return(target/all_sum)
}
probpre <-function(param, k, x, t){
  all_sum <- 0
  tmp <- 0
  for(m in 1:length(param$category_prm)){
    tmp <- tmp + 1.7 * logit(param, param$category_prm[m], x, t)
    if(m==k) target <- exp(tmp)
    all_sum <- all_sum + exp(tmp)
  }
  return(target/all_sum)
}
gen_category_param <- function(N, K){
  category <- matrix(0, nrow=N, ncol=(K-1))
  for (k in 1:(K-1)){
    category[,k] <- rnorm(N, 0, 1)
  }
  for (i in 1:N){
    category[i,] = sort(category[i,])
    category[i,] = category[i,] - mean(category[i,])
  } 
  return(category)
}

gen_category_param_previous <- function(K){
  category <- matrix(0, nrow=1, ncol=(K-1))
  for (k in 1:(k-1)){
    category[k] <- rnorm(1, 0, 1)
  }
  category = sort(category)
  category <- category - mean(category)
  return(category)
}

read_data <- function(setting, filename){
  data <- read.table(filename, header=TRUE,sep=",")
  colnames(data) <- c("ExamineeID", "Score","RaterID","TimeID")
  data_irt=list(
    K=setting$K, J=setting$n_person, R=setting$n_rater, T=setting$n_time, N=nrow(data), 
    ExamineeID=data$ExamineeID, RaterID=data$RaterID, TimeID=data$TimeID, X=data$Score)
  return(data_irt)
}

get_alpha_estimates_with_restriction <- function(d){
  return(append(1.0/prod(d), d))
}

get_estimates_with_mean_restriction <- function(d){
  return(append(-1*sum(d), d))
}

log_normal_generator <- function(n, mean, sd) {
  sdlog <- sqrt(log((sd/mean)^2 + 1))
  meanlog <- log(mean) - (sdlog^2) / 2
  rlnorm(n, meanlog = meanlog, sdlog = sdlog)
}

calc_rmse <- function(a, b){
  return(sqrt(mean((a - b)^2)))
}

calc_bias <- function(a, b){
  return(mean(a - b))
}

calc_error <- function(a, b){
  return(list(RMSE = calc_rmse(a, b), BIAS = calc_bias(a, b)))
}

generate_constrained_alpha <- function(N){
  const_alpha <- log_normal_generator(N, mean=1, sd=0.4)
  const_alpha[2:N] = exp(log(const_alpha[2:N]) - mean(log(const_alpha[2:N])))
  const_alpha[1] = 1.0
  return(const_alpha)
}

generate_beta_rt <- function(R,T){
  const_beta_rt <- matrix(0, nrow=R, ncol=T)
  for (r in 1:R){
    const_beta_rt[r,1] <- rnorm(1);
    for (t in 2:T){
      const_beta_rt[r,t] <- rnorm(1, const_beta_rt[r,t-1], 0.1);
    }
  }
  return(const_beta_rt)
}

generate_beta_rt_self <- function(R,T){
  const_beta_rt <- matrix(0, nrow=R, ncol=T);
  r2 = trunc(R/3);
  r3 = r2 * 2;
  t2 = trunc(T/2)
  for (r in 1:r2){
    const_beta_rt[r,1] <- rnorm(1, 0, 1);
    const_beta_rt[r,2:t2] <- const_beta_rt[r,1];
    const_beta_rt[r,(t2+1)] <-const_beta_rt[r,1]+rnorm(1,0,0.2);
    const_beta_rt[r,(t2+1):T] <- const_beta_rt[r,(t2+1)];
  }
  for (r in (r2+1):r3){
    const_beta_rt[r,1] <- rnorm(1, 0, 1);
    for(t in 2:T){
      const_beta_rt[r,t] <- rnorm(1, 0, 0.01) + const_beta_rt[r,1];
    }
  }
  for (r in (r3+1):R){
    const_beta_rt[r,1] <- rnorm(1, 0, 1);
    for(t in 2:T){
    const_beta_rt[r,t] <- 1.1 * const_beta_rt[r,(t-1)];
    }
  }
  return(const_beta_rt)
}

generate_data <- function(param, setting){
  N <- setting$n_item * setting$n_person * setting$n_rater
  U = matrix(0, nrow=N, ncol=5)
  tsub = setting$n_person/setting$n_time
  if(setting$n_person%%setting$n_time != 0){
    tsub = trunc(setting$n_person/setting$n_time)
  }
  tsub2 = 0
  tsub3 = 0
  row_idx = 1
  for(t in 1:setting$n_time){
    tsub2 = tsub * t - (tsub - 1)
    tsub3 = tsub * t
    if(tsub3 > setting$n_person){
      tsub3 = setting$n_person
    }
    if(t == setting$n_time && tsub3 < setting$n_person){
      tsub3 = setting$n_person
    }
      for (j in tsub2:tsub3){
        for (i in 1:setting$n_item){
          for (r in 1:setting$n_rater){
              prob <- c(1:setting$K)
              for (k in 1:(setting$K)){
                prob[k] = prob(get_prm_list(param, i, r, t), k, param$theta[j])
              }
              score = grep(1, rmultinom(1, 1, prob))
              U[row_idx,] <- c(j, i, r, t, score)
              row_idx = row_idx + 1
          }
        }
      
    }
  }
  colnames(U) <- c("ExamineeID", "ItemID","RaterID","TimeID","Score")
  data <- data.frame(U)
  data=list(K=setting$K, J=setting$n_person, I=setting$n_item, R=setting$n_rater, T=setting$n_time, N=nrow(data), 
            ItemID=data$ItemID, ExamineeID=data$ExamineeID, RaterID=data$RaterID, TimeID=data$TimeID, X=data$Score)
  return(data)
}

get_theta_se_stat <- function(fit){
  SE <- summary(fit, par="theta")$summary[,"sd"]
  return(list(mean_se = mean(SE), sd_se = sd(SE)))
}

draw_icc <- function(param, K, def){
  x<- seq(-def$xlim, def$xlim, length=(10))
  curve(prob(param, 1, x), 
        from=-def$xlim, to=def$xlim, ylim=c(0,def$ylim), xlab=def$xlab, ylab=def$ylab, main=def$title, 
        col=def$color[1],lty=def$style[1], cex.lab=def$lcex, lwd=def$llwd, cex.axis=def$axcx, cex.main=def$maincx, yaxt= "n")
  par(new=T)
  plot(x, prob(param, 1, x), xlim=c(-def$xlim, def$xlim), ylim=c(0,def$ylim), xlab="", ylab="",
       col=def$color[1],lty=def$style[1], cex.lab=def$lcex, lwd=def$llwd, axes=FALSE, pch=def$pchs[1], cex=1.5,  yaxt= "n")
  for(k in 2:K){
    par(new=T)
    curve(prob(param, k, x), 
          from=-def$xlim, to=def$xlim, ylim=c(0,def$ylim), xlab="", ylab="",  
          col=def$color[k],lty=def$style[k], cex.lab=def$lcex, lwd=def$llwd, axes=FALSE,  yaxt= "n" )
    par(new=T)
    plot(x, prob(param, k, x), xlim=c(-def$xlim, def$xlim), ylim=c(0,def$ylim), xlab="", ylab="",
         col=def$color[k],lty=def$style[k], cex.lab=def$lcex, lwd=def$llwd, axes=FALSE, pch=def$pchs[k], cex=1.5,  yaxt= "n")
  }
  axis(2, cex.axis = def$caxcx)
}
