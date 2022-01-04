data{
  int <lower=0> J;//n_examinee
  int <lower=0> R;//n_rater
  int <lower=2> K;//n_score
  int <lower=0> T;//n_time
  int <lower=0> N;//n_samples
  int <lower=1, upper=J> ExamineeID [N];
  int <lower=1, upper=R> RaterID [N];
  int <lower=1, upper=K> X [N];//Score
  int <lower=1, upper=T> TimeID[N];
}
transformed data{
  vector[K] c = cumulative_sum(rep_vector(1, K)) - 1;
}

parameters {
  vector[J] theta;
  real<lower=0> alpha_r [R-1];
  matrix[R,T] beta_rt;
  vector[K-2] beta_rk;
  real<lower=0> sigma_beta_rt;
}
transformed parameters{
  vector[K-1] category_est;
  vector[K] category_prm;
  real<lower=0> trans_alpha_r[R];
  trans_alpha_r[1] = 1.0 / prod(alpha_r);
  trans_alpha_r[2:R] = alpha_r;
  category_est[1:(K-2)] = beta_rk;
  category_est[K-1] = -1*sum(beta_rk);
  category_prm = cumulative_sum(append_row(0, category_est));
}
model{
  theta ~ normal(0, 1);
  trans_alpha_r ~ lognormal(0.0, 0.4);
  sigma_beta_rt ~ lognormal(-3, 1);
  for (r in 1:R){
    beta_rt[r,1] ~ normal(0, 1);
    for (t in 2:T){
      beta_rt[r,t] ~ normal(beta_rt[r,t-1], sigma_beta_rt);
    }
  }
  category_est  ~ normal(0, 1);
  
  for (n in 1:N){
    target += 1/log(N) * categorical_logit_log(X[n],1.7*trans_alpha_r[RaterID[n]]*(c*(theta[ExamineeID[n]]-beta_rt[RaterID[n],TimeID[n]])-category_prm));
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] = categorical_logit_log(X[n], 1.7*trans_alpha_r[RaterID[n]]*(c*(theta[ExamineeID[n]]-beta_rt[RaterID[n],TimeID[n]])-category_prm));
  }
}
