data {
  int N;       // i=1:N observations
  int M;       // j=1:M levels
  vector[N] y;
  array[N] int group;
}
parameters {
  vector[M] b;
  real<lower=0> sigma;
}
model {
  for(j in 1:M){
    b[j] ~ normal(25,10); 
  }
  // or short: b ~ normal(25,10);
  sigma ~ exponential(0.1);
  for(i in 1:N){
    y[i] ~ normal(b[group[i]], sigma);
  }
  // or short: y ~ normal(b[group], sigma);
}
generated quantities{
  vector[N] log_lik;
  for(i in 1:N){
    log_lik[i] = normal_lpdf(y[i] | b[group[i]], sigma);
  }
}
