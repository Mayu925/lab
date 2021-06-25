get_theta_se_stat <- function(fit){
  SE <- summary(fit, par="theta")$summary[,"sd"]
  return(list(mean_se = mean(SE), sd_se = sd(SE)))
}

generate_data <- function(param, setting){
  N <- setting$n_item * setting$n_person * setting$n_rater * setting$n_time
  U = matrix(0, nrow=N, ncol=5)
  row_idx = 1
  for (j in 1:setting$n_person){
    for (i in 1:setting$n_item){
      for (r in 1:setting$n_rater){
        for (t in 1:setting$n_time){
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

get_eap_each <- function(data, param){
  dat <- data
  Xh <- seq(-3, 3, length=150)
  AXh <- dnorm(Xh, 0, 1)
  for(n in 1:data$N){
    Z1 = 0.0;
    Z2 = 0.0;
    for(h in 1:length(Xh)){
      if(is.na(data$X[n]) == FALSE){
        p <- prob(get_prm_list(param, data$ItemID[n], data$RaterID[n], data$TimeID[n]),
                data$X[n], Xh[h])              
        Z1 = Z1 + exp(p)*AXh[h]*Xh[h];
        Z2 = Z2 + exp(p)*AXh[h];
      }
    }
    dat$X[n] = (Z1 / Z2);
  }
  return(dat)
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
