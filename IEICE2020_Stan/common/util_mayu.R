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
