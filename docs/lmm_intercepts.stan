data {
  int N;       // i=1:N observations
  int M;       // j=1:M levels
  vector[N] y;
  array[N] int group;
}
parameters {
  real mu_b;
  real<lower=0> sd_b;
  vector[M] b;
  real<lower=0> sigma;
}
model {
  for(j in 1:M){
    b[j] ~ normal(mu_b,sd_b); 
  }
  // or short: b ~ normal(mu_b,sd_b);
  mu_b ~ normal(25,10);
  sd_b ~ exponential(0.1);
  sigma ~ exponential(0.1);
  for(i in 1:N){
    y[i] ~ normal(b[group[i]], sigma);
  }
  // or short: y ~ normal(b[group], sigma);
}
