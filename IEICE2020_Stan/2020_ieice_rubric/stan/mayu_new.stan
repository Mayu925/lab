data{
  int <lower=0> J;//n_examinee
  int <lower=0> I;//n_item
  int <lower=0> R;//n_rater
  int <lower=2> K;//n_score
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
  real<lower=0> alpha_r [R];
  vector[K-2] beta_rk[R];
  matrix[R,T] alpha_rt;

}
transformed parameters{
  vector[K-1] category_est[R];
  vector[K] category_prm[R];

  for(p in 1:R){
    category_est[p, 1:(K-2)] = beta_rk[p];
    category_est[p, K-1] = -1*sum(beta_rk[p]);
      category_prm[p] = cumulative_sum(append_row(0, category_est[p]));
  }
  
}

model{
  alpha_r ~ lognormal(0.0, 1.0);
  theta ~ normal(0, 1);

  
  for (r in 1:R){
    alpha_rt[r,1] ~ normal(0, 1);
    for (t in 2:T){
      alpha_rt[r,t] ~ normal(alpha_rt[r,t-1], 1);
    }
  }
  for (p in 1:R) category_est [p,] ~ normal(0, 1);
  
  for (n in 1:N){
    X[n] ~ categorical_logit(alpha_r[RaterID[n]]*(c*(theta[ExamineeID[n]])-alpha_rt[RaterID[n],TimeID[n]]-category_prm[RaterID[n]]));
  }
}


generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] = categorical_logit_log(X[n], alpha_r[RaterID[n]]*(c*(theta[ExamineeID[n]])-alpha_rt[RaterID[n],TimeID[n]]-category_prm[RaterID[n]]));
  }
}
