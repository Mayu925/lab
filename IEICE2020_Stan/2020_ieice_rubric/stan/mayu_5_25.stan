data{
  int <lower=0> J;//n_examinee
  int <lower=0> I;//n_item
  int <lower=0> R;//n_rater
  int <lower=0> C;//n_rubric
  int <lower=2> K;
  int <lower=0> T;//n_time
  int <lower=0> N;//n_samples
  int <lower=1, upper=J> ExamineeID [N];
  int <lower=1, upper=I> ItemID [N];
  int <lower=1, upper=R> RaterID [N];
  //int <lower=1, upper=C> RubricID [N];
  int <lower=1, upper=K> X [N];//Score
  int <lower=1, upper=T> TimeID[N];
}
transformed data{
  vector[K] c = cumulative_sum(rep_vector(1, K)) - 1;
}
parameters {
  vector[J] theta;
  real<lower=0> alpha_i [I-1];
  real<lower=0> alpha_r [R];
  real<lower=0> alpha_c [C-1];
  vector[I-1] beta_i;
  vector[R] beta_r;
  vector[C-1] beta_c;
  vector[K-2] beta_ck [C];
  vector[K] beta_rk [R];
  vector[K] severity[R, T];
  real<lower=0> tau_r [R-1];
  vector[T] alpha_rt[R];
  real<lower=0> pai_1r[R];
}
transformed parameters{
  real<lower=0> trans_alpha_i[I];
  real<lower=0> trans_alpha_c[C];
  vector[I] trans_beta_i;
  vector[C] trans_beta_c;
  real<lower=0> trans_tau_r[R];
  vector[K-1] category_est[C];
  vector[K] category_prm[C, R];
  //vector[K] severity[R, T];
  trans_alpha_i[1] = 1.0 / prod(alpha_i);
  trans_alpha_c[1] = 1.0 / prod(alpha_c);
  trans_beta_i[1] = -1*sum(beta_i);
  trans_beta_c[1] = -1*sum(beta_c);
  trans_tau_r[1] = 1.0 / prod(tau_r);
  trans_alpha_i[2:I] = alpha_i;
  trans_alpha_c[2:C] = alpha_c;
  trans_beta_i[2:I] = beta_i;  
  trans_beta_c[2:C] = beta_c;  
  trans_tau_r[2:R] = tau_r;  
  for(p in 1:C){
    category_est[p, 1:(K-2)] = beta_ck [p];
    category_est[p, K-1] = -1*sum(beta_ck [p]);  
    for(r in 1:R){
      category_prm[p, r] = cumulative_sum(append_row(0, trans_tau_r[r]*category_est[p]));//τ_r*d_cm
    }
    
  }
  /*
  for(t in 1:T){
    for(r in 1:R){
      severity[r,t] = cumulative_sum(append_row(0, pai_1r[r]*alpha_rt[t])); 
    }
  }*/
}
model{
  trans_alpha_i ~ lognormal(0.0, 1.0);
  alpha_r ~ lognormal(0.0, 1.0);
  trans_alpha_c ~ lognormal(0.0, 1.0);
  trans_beta_i ~ normal(0, 1);
  beta_r ~ normal(0, 1);
  trans_beta_c ~ normal(0, 1.0);
  trans_tau_r ~ lognormal(0.0, 0.5);
  theta ~ normal(0, 1);
  for (r in 1:R) severity[r,1] ~ normal(0, 1);
  for (p in 1:C) category_est [p,] ~ normal(0, 1);
  for (n in 1:N){
    X[n] ~ categorical_logit(alpha_r[RaterID[n]]*(theta[ExamineeID[n]]-severity[RaterID[n],TimeID[n]]-beta_rk[RaterID[n]]));
  }
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] = categorical_logit_log(X[n],alpha_r[RaterID[n]]*(theta[ExamineeID[n]]-severity[RaterID[n],TimeID[n]]-beta_rk[RaterID[n]]));
  }
}
