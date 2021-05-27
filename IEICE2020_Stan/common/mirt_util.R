get_theta_se_stat <- function(fit){
  SE <- summary(fit, par="theta")$summary[,"sd"]
  return(list(mean_se = mean(SE), sd_se = sd(SE)))
}

generate_data <- function(param, setting){
  N <- setting$n_item * setting$n_person * setting$n_rater * setting$n_rubric
  U = matrix(0, nrow=N, ncol=5)
  row_idx = 1
  for (j in 1:setting$n_person){
    for (i in 1:setting$n_item){
      for (r in 1:setting$n_rater){
        for (c in 1:setting$n_rubric){
          prob <- c(1:setting$K)
          for (k in 1:(setting$K)){
            prob[k] = prob(get_prm_list(param, i, r, c), k, param$theta[j,])
          }
          score = grep(1, rmultinom(1, 1, prob))
          U[row_idx,] <- c(j, i, r, c, score)
          row_idx = row_idx + 1
        }
      }
    }
  }
  colnames(U) <- c("ExamineeID", "ItemID","RaterID","RubricID","Score")
  data <- data.frame(U)
  data=list(K=setting$K, J=setting$n_person, I=setting$n_item, R=setting$n_rater, C=setting$n_rubric, N=nrow(data), D=setting$n_dim,
            ItemID=data$ItemID, ExamineeID=data$ExamineeID, RaterID=data$RaterID, RubricID=data$RubricID, X=data$Score)
  return(data)
}

get_likelihood <- function(data, param){
  X <- to_wide_format(data)
  LL <- 0
  for(j in 1:data$J){
    for(i in 1:data$I){
      for(r in 1:data$R){
        for(c in 1:data$C){
          k <- X[j, i, r, c]
          if(is.na(k) == FALSE){
            p <- prob(get_prm_list(param, i, r, c), k, param$theta[j,])              
            LL = LL + log(p) 
          }
        }
      }
    }
  }
  return (LL)       
}
