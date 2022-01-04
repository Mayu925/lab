data{
  int <lower=0> J;//n_examinee
  int <lower=0> R;//n_rater
  int <lower=2> K;
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
  vector[K-2] beta_rk[R];
  vector[R] pai_1r;
  vector[R-1] pai_0r;
  real<lower=0> alpha_r [R-1];
}
transformed parameters{
  vector[K-1] category_est[R];
  vector[K] category_prm[R];
  real<lower=0> trans_alpha_r[R];
  vector[R] trans_pai_0r;
  trans_pai_0r[1:(R-1)] = pai_0r;
  trans_pai_0r[R] = -1*sum(pai_0r); 
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
  pai_1r ~ normal(0, 1);
  trans_alpha_r ~ lognormal(0.0, 0.4);
  pai_0r ~ normal(0, 1);
  for (p in 1:R) category_est [p,] ~ normal(0, 1);
  
  pai_1r ~ normal(0, 1);
  trans_pai_0r ~ normal(0, 1);
  for (n in 1:N){
    X[n] ~ categorical_logit(1.7*trans_alpha_r[RaterID[n]]*(c*(theta[ExamineeID[n]]-trans_pai_0r[RaterID[n]]-pai_1r[RaterID[n]]*(TimeID[n]/T))-category_prm[RaterID[n]]));
  }
}


generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] = categorical_logit_log(X[n], 1.7*trans_alpha_r[RaterID[n]]*(c*(theta[ExamineeID[n]]-trans_pai_0r[RaterID[n]]-pai_1r[RaterID[n]]*(TimeID[n]/T))-category_prm[RaterID[n]]));
  }
}
