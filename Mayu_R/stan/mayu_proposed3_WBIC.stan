data{
  int <lower=0> J;//n_examinee
  int <lower=0> I;//n_item
  int <lower=0> R;//n_rater
  int <lower=2> K;
  int <lower=0> T;//n_time
  int <lower=0> N;//n_samples
  int <lower=1, upper=J> ExamineeID [N];
  int <lower=1, upper=I> ItemID [N];
  int <lower=1, upper=R> RaterID [N];
  int <lower=1, upper=K> X [N];//Score
  int <lower=1, upper=T> TimeID[N];
}
transformed data{
  vector[K] c = cumulative_sum(rep_vector(1, K)) - 1;
}

parameters {
  vector[J] theta;
  vector[K-2] beta_rk[R];
  matrix[R,T] beta_rt;
  real<lower=0> pai_1r [R-1];
  vector[R] pai_0r;
  real<lower=0> alpha_r [R-1];
  real<lower=0> sigma_beta_rt;
}
transformed parameters{
  vector[K-1] category_est[R];
  vector[K] category_prm[R];
  real<lower=0> trans_pai_1r[R];
  real<lower=0> trans_alpha_r[R];
  trans_pai_1r[1] = 1.0 / prod(pai_1r);
  trans_pai_1r[2:R] = pai_1r;
  trans_alpha_r[1] = 1.0 / prod(alpha_r);
  trans_alpha_r[2:R] = alpha_r;
  for(r in 1:R){
    category_est[r, 1:(K-2)] = beta_rk[r];
    category_est[r, K-1] = -1*sum(beta_rk[r]);
    category_prm[r] = cumulative_sum(append_row(0, category_est[r]));
  }
}


model{
  theta ~ normal(0, 1);
  trans_pai_1r ~ lognormal(0.0, 1.0);
  trans_alpha_r ~ lognormal(0.0, 0.4);
  pai_0r ~ normal(0, 1);
  sigma_beta_rt ~ lognormal(-3, 1);
  for (r in 1:R){
    beta_rt[r,] ~ normal(0, 1);
  }
   for (p in 1:R) category_est [p,] ~ normal(0, 1);
  
  for (n in 1:N){
    target += 1/log(N) * categorical_logit_log(X[n],1.7*trans_alpha_r[RaterID[n]]*(c*(theta[ExamineeID[n]]-pai_0r[RaterID[n]]-trans_pai_1r[RaterID[n]]*beta_rt[RaterID[n],TimeID[n]])-category_prm[RaterID[n]]));
  }
}


generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] = categorical_logit_log(X[n], 1.7*trans_alpha_r[RaterID[n]]*(c*(theta[ExamineeID[n]]-pai_0r[RaterID[n]]-trans_pai_1r[RaterID[n]]*beta_rt[RaterID[n],TimeID[n]])-category_prm[RaterID[n]]));
  }
}
