get_result_statistics_common <- function(fit1, data1, setting){
  est_param <- get_estimates(fit1, setting)
  Rhat <- get_Rhat_stat(fit1)
  SE_stat <- get_theta_se_stat(fit1)
  D <- c("WAIC", waic(extract_log_lik(fit1))$estimates[3,1], 
         "ML", get_ml(fit1),
         "AIC", aic(data1, est_param),
         "BIC", bic(data1, est_param),
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
    prm = category_prm[1:K-1]
    prm = append(prm, -1*sum(prm))
    mat = t(data.frame(prm))
  return(mat)
}
convert_alpha_rt <- function(alpha_rt, T, R){
  for(r in 1:R){
    alrt = alpha_rt[((r-1)*T+1):((r-1)*T+T)]
    if(r == 1){
      mat = t(data.frame(alrt))
    }else{
      mat = rbind(mat,t(data.frame(alrt)))
    }
  }
  return(mat)
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
  for (k in 1:(K-1)){
    category <- rnorm(K, 0, 1)
  }
  return(category)
}


read_data <- function(setting, filename){
  data <- read.table(filename, header=TRUE,sep=",")
  colnames(data) <- c("ExamineeID", "ItemID", "RaterID", "Score", "TimeID")
  data_irt=list(
    K=setting$K, J=setting$n_person, I=setting$n_item, R=setting$n_rater, T=setting$n_time, N=nrow(data), 
    ItemID=data$ItemID, ExamineeID=data$ExamineeID, RaterID=data$RaterID, TimeID=data$TimeID, X=data$Score)
  return(data_irt)
}

get_alpha_estimates_with_restriction <- function(d){
  return(append(1.0/prod(d), d))
}

get_estimates_with_mean_restriction <- function(d){
  return(append(-1*sum(d), d))
}

aic <- function(data, param){
  return(-2*get_likelihood(data, param) + 2*(get_param_size(data)));
}

bic <- function(data, param){
  return(-2*get_likelihood(data, param) + (get_param_size(data)) *log(data$J));
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

to_wide_format <- function(data){
  X = array(NA, dim = c(data$J, data$I, data$R, data$T))
  for (n in 1:length(data$X)){
    X[data$ExamineeID[n], data$ItemID[n], data$RaterID[n],data$TimeID[n]] = data$X[n];
  }
  return(X)
}

to_long_format <- function(data, X){
  for (n in 1:length(data$X)){
    data$X[n] <- X[data$ExamineeID[n], data$ItemID[n], data$RaterID[n],data$TimeID[n]];
  }
  return(data)
}

generate_constrained_alpha <- function(N){
  const_alpha <- log_normal_generator(N, mean=1, sd=0.4)
  const_alpha[2:N] = exp(log(const_alpha[2:N]) - mean(log(const_alpha[2:N])))
  const_alpha[1] = 1.0
  return(const_alpha)
}

generate_alpha_rt <- function(R,T){
  for (r in 1:R){
    const_alpha_rt[r,1] <- rnorm(1);
    for (t in 2:T){
      const_alpha_rt[r,t] <- rnorm(1, alpha_rt[r,t-1], 1);
    }
  }
  return(const_alpha_rt)
}

get_ml <- function(fit){
  ms <- rstan::extract(fit)
  theta <- apply(ms$log_lik, 1, sum)
  return(-2*(1/sum(1/theta)*length(theta)))
}

generate_data <- function(param, setting){
  N <- setting$n_item * setting$n_person * setting$n_rater
  U = matrix(0, nrow=N, ncol=5)
  row_idx = 1
  for (j in 1:setting$n_person){
    for (i in 1:setting$n_item){
      for (r in 1:setting$n_rater){
        for (t in 1:setting$n_time){
          if(j==t){
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
  }
  colnames(U) <- c("ExamineeID", "ItemID","RaterID","TimeID","Score")
  data <- data.frame(U)
  data=list(K=setting$K, J=setting$n_person, I=setting$n_item, R=setting$n_rater, T=setting$n_time, N=nrow(data), 
            ItemID=data$ItemID, ExamineeID=data$ExamineeID, RaterID=data$RaterID, TimeID=data$TimeID, X=data$Score)
  return(data)
}

get_eap <- function(X, data, param){
  Xh <- seq(-3, 3, length=50)
  AXh <- dnorm(Xh, 0, 1)
  theta_eap <- c(1:data$J)
  names(theta_eap) <- c(1:data$J)
  for(j in 1:data$J){
    Z1 = 0.0;
    Z2 = 0.0;
    for(h in 1:length(Xh)){
      LL = 0.0;
      for(i in 1:data$I){
        for(r in 1:data$R){
          for(t in 1:data$T){
            k <- X[j, i, r, t]
            if(is.na(k) == FALSE){
              p <- prob(get_prm_list(param, i, r, t), k, Xh[h])              
              LL = LL + log(p);
            }
          }
        }
      }
      Z1 = Z1 + exp(LL)*AXh[h]*Xh[h];
      Z2 = Z2 + exp(LL)*AXh[h];
    }
    theta_eap[j] = (Z1 / Z2);
  }  
  return(theta_eap)
}

get_likelihood <- function(data, param){
  X <- to_wide_format(data)
  LL <- 0
  for(j in 1:data$J){
    for(i in 1:data$I){
      for(r in 1:data$R){
        for(t in 1:data$T){
          k <- X[j, i, r, t]
          if(is.na(k) == FALSE){
            p <- prob(get_prm_list(param, i, r, t), k, param$theta[j])              
            LL = LL + log(p) 
          }
        }
      }
    }
  }
  return (LL)       
}

get_theta_se_stat <- function(fit){
  SE <- summary(fit, par="theta")$summary[,"sd"]
  return(list(mean_se = mean(SE), sd_se = sd(SE)))
}
