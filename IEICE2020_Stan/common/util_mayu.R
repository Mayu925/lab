read_data <- function(setting, filename){
  data <- read.table(filename, header=TRUE,sep=",")
  colnames(data) <- c("ExamineeID", "ItemID", "RaterID", "TimeID", "Score")
  data_irt=list(
    K=setting$K, J=setting$n_person, I=setting$n_item, R=setting$n_rater, T=setting$n_rubric, N=nrow(data), 
    ItemID=data$ItemID, ExamineeID=data$ExamineeID, RaterID=data$RaterID, TimeID=data$TimeID, X=data$Score)
  return(data_irt)
}

get_alpha_estimates_with_restriction <- function(d){
  return(append(1.0/prod(d), d))
}

get_estimates_with_mean_restriction <- function(d){
  return(append(-1*sum(d), d))
}

standard_error <- function(param, x){
  return(sqrt(1/fisher_information(param, x)))
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
  X = array(NA, dim = c(data$J, data$I, data$R, data$C))
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

get_average_rubric <- function(data){
  X <- to_wide_format(data)
  TH <- array(-1, dim = c(data$J, data$C))
  for(j in 1:data$J){
    for(c in 1:data$C){  
      TH[j,c] <- mean(X[j, , , c])
    }
  }
  return(TH)
}

get_average_rater <- function(data){
  X <- to_wide_format(data)
  TH <- array(-1, dim = c(data$J, data$R))
  for(j in 1:data$J){
    for(r in 1:data$R){  
      TH[j,r] <- mean(X[j, , r, ])
    }
  }
  return(TH)
}

get_average_item <- function(data){
  X <- to_wide_format(data)
  TH <- array(-1, dim = c(data$J, data$I))
  for(j in 1:data$J){
    for(i in 1:data$I){
      TH[j,i] <- mean(X[j, i, , ])
    }
  }
  return(TH)
}

get_rubric_FI <- function(param, i, r, set_C, x){
  inf <- 0
  for(c in set_C)
    inf <- inf + fisher_information(get_prm_list(param, i=i, r=r, c=c), x)
  return(inf)        
}

get_rater_FI <- function(param, i, set_R, c, x){
  inf <- 0
  for(r in set_R)
    inf <- inf + fisher_information(get_prm_list(param, i=i, r=r, c=c), x)
  return(inf)        
}

generate_constrained_alpha <- function(N){
  const_alpha <- log_normal_generator(N, mean=1, sd=0.4)
  const_alpha[2:N] = exp(log(const_alpha[2:N]) - mean(log(const_alpha[2:N])))
  const_alpha[1] = 1.0
  return(const_alpha)
}

get_eap_rubric <- function(data, param){
  TH <- array(-1, dim = c(data$J, data$C))
  for(c in 1:data$C){
    X <- to_wide_format(data)
    X[,,,c(1:data$C)[-c]] <- NA
    TH[, c] = get_eap(X, data, param);
  }
  return(TH)
}

get_eap_rater <- function(data, param){
  TH <- array(-1, dim = c(data$J, data$R))
  for(r in 1:data$R){
    X <- to_wide_format(data)
    X[,,c(1:data$R)[-r],] <- NA
    TH[, r] = get_eap(X, data, param);
  }  
  return(TH)
}

get_eap_item <- function(data, param){
  TH <- array(-1, dim = c(data$J, data$I))
  for(i in 1:data$I){
    X <- to_wide_format(data)
    X[,c(1:data$I)[-i],,] <- NA
    TH[, i] = get_eap(X, data, param);
  }  
  return(TH)
}

get_ml <- function(fit){
  ms <- rstan::extract(fit)
  theta <- apply(ms$log_lik, 1, sum)
  return(-2*(1/sum(1/theta)*length(theta)))
}

get_generalizability <- function(data){
  library(lme4) 
  dat <- data
  dat$ItemID <- factor(dat$ItemID)
  dat$ExamineeID <- factor(dat$ExamineeID)
  dat$RubricID <- factor(dat$RubricID)
  dat$RaterID <- factor(dat$RaterID)
  mod1 <- lmer(X ~ 1 + (1|ItemID) + (1|ExamineeID) + (1|RubricID)  + (1|RaterID) 
               + (1|ExamineeID:ItemID) + (1|ExamineeID:RubricID) + (1|ExamineeID:RaterID)
               + (1|ItemID:RubricID) + (1|ItemID:RaterID) + (1|RubricID:RaterID), data=dat)
#  mod1 <- lmer(X ~ 1 + (1|ItemID) + (1|ExamineeID) + (1|RubricID)  + (1|RaterID) 
#               + (1|ExamineeID:ItemID) + (1|ExamineeID:RubricID) + (1|ExamineeID:RaterID), data=dat)

  # VarCorrでrandom effectsの出力を取り出す
  varcomp <- VarCorr(mod1)  
  output <- rep(NA, length(varcomp)+1)
  for(i in 1:length(varcomp)){
    output[i] <- varcomp[[i]][,1]
  }
  names(output) <- c(names(varcomp), "Residual")
  output <- sort(output, na.last=T, decreasing=T)
  output["Residual"] <- attr(varcomp, "sc")[[1]]^2
  print(round(output/sum(output),3)*100)
  
  p  <- output[["ExamineeID"]]
  pr <- output[["ExamineeID:RaterID"]]
  pi <- output[["ExamineeID:ItemID"]]
  pc <- output[["ExamineeID:RubricID"]]
  pri  <- output[["Residual"]]
  denom <- p + pr / data$R + pi / data$I + pc / data$C + pri / (data$R * data$I * data$C)
  G  <- p / denom
  return(G)
}
