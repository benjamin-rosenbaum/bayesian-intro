data {
  int N;       // i=1:N observations
  int M;       // j=1:M levels
  vector[N] y;
  vector[N] x; // continuous predictor
  array[N] int group; // categorical predictor
}
parameters {
  real b[M];   // M intercepts
  real c;      // joint slope
  real<lower=0> sigma;
}
model {
  for(j in 1:M){
    b[j] ~ normal(25,10); 
  }
  // or short: b ~ normal(25,10);
  c ~ normal(5, 10);
  sigma ~ exponential(0.1);
  for(i in 1:N){
    y[i] ~ normal(b[group[i]]+c*x[i], sigma);
  }
  // or short: y ~ normal(b[group]+c*x, sigma);
}
