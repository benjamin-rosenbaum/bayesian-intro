data {
  int<lower=0> N;
  vector[N] dose;
  array[N] int total;
  array[N] int y;
}
parameters {
  real a;
  real b;
}
model {
  real mu;
  a ~ normal(0,1);
  b ~ normal(1,1);
  for(i in 1:N){
    mu = inv_logit(a+b*dose[i]);
    y[i] ~ binomial(total[i], mu);
  }
  // short: y ~ binomial_logit(total, dose, a, b);
}
