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
  vector[K-2] beta_rk;
  matrix[R,T] alpha_rt;
  real<lower=0> pai_1r [R];
  vector[R] pai_0r;

}
transformed parameters{
  vector[K-1] category_est;
  vector[K] category_prm;

    category_est[1:(K-2)] = beta_rk;
    category_est[K-1] = -1*sum(beta_rk);
      category_prm = cumulative_sum(append_row(0, category_est));
}

model{
  theta ~ normal(0, 1);
  pai_1r ~ lognormal(0.0, 1.0);
  pai_0r ~ normal(0, 1);
  
  for (r in 1:R){
    alpha_rt[r,1] ~ normal(0, 1);
    for (t in 2:T){
      alpha_rt[r,t] ~ normal(0, 1);
    }
  }
  category_est  ~ normal(0, 1);
  
  for (n in 1:N){
    X[n] ~ categorical_logit(c*(theta[ExamineeID[n]])-pai_0r[RaterID[n]]-pai_1r[RaterID[n]]*alpha_rt[RaterID[n],TimeID[n]]-category_prm);
  }
}


generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] = categorical_logit_log(X[n], c*(theta[ExamineeID[n]])-pai_0r[RaterID[n]]-pai_1r[RaterID[n]]*alpha_rt[RaterID[n],TimeID[n]]-category_prm);
  }
}
