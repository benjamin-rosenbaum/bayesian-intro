rm(list=ls())
library(rstan)
library(coda)

setwd("~/Desktop/teaching Bayes")


set.seed(123) # initiate random number generator for reproducability
n=50

a=1
b=1
sigma=1

x = runif(n=n, min=-1, max=1)
y = rnorm(n=n, mean=a+b*x, sd=sigma)

df = data.frame(x=x,
                y=y)

plot(df)

stan_code_1 = '
data {
  int n;
  vector[n] x;
  vector[n] y;
}
parameters {
  vector[2] b;  
  real<lower=0> sigma;  // standard deviation
}
model {
  // priors
  b[1] ~ normal(0, 10);
  b[2] ~ normal(0, 10);
  sigma ~ uniform(0, 100);
  // likelihood
  y ~ normal(b[1]+b[2]*x, sigma);
}
'

stan_code_2 = '
data {
  int n;
  vector[n] x;
  vector[n] y;
}
parameters {
  vector[2] b;  
  real<lower=0> sigma;  // standard deviation
}
model {
  // priors
  b[1] ~ normal(0, 10);
  b[2] ~ normal(0, 1.2);
  sigma ~ uniform(0, 100);
  // likelihood
  y ~ normal(b[1]+b[2]*x, sigma);
}
'

stan_code_3 = '
data {
  int n;
  vector[n] x;
  vector[n] y;
}
parameters {
  vector[2] b;  
  real<lower=0> sigma;  // standard deviation
}
model {
  // priors
  b[1] ~ normal(0, 10);
  b[2] ~ normal(0, 0.5);
  sigma ~ uniform(0, 100);
  // likelihood
  y ~ normal(b[1]+b[2]*x, sigma);
}
'

stan_code_4 = '
data {
  int n;
  vector[n] x;
  vector[n] y;
}
parameters {
  vector[2] b;  
  real<lower=0> sigma;  // standard deviation
}
model {
  // priors
  b[1] ~ normal(0, 10);
  b[2] ~ normal(0, 0.15);
  sigma ~ uniform(0, 100);
  // likelihood
  y ~ normal(b[1]+b[2]*x, sigma);
}
'

data = list(n=n, 
            x=df$x, 
            y=df$y)

rstan_options(auto_write = TRUE)
options(mc.cores = 4) 

stan_model_1 = stan_model(model_code=stan_code_1)
stan_model_2 = stan_model(model_code=stan_code_2)
stan_model_3 = stan_model(model_code=stan_code_3)
stan_model_4 = stan_model(model_code=stan_code_4)

fit_1  = sampling(stan_model_1,
                  data=data,
                  chains=4,
                  iter=30000,
                  warmup=1000
)
fit_2  = sampling(stan_model_2,
                  data=data,
                  chains=4,
                  iter=30000,
                  warmup=1000
)
fit_3  = sampling(stan_model_3,
                  data=data,
                  chains=4,
                  iter=30000,
                  warmup=1000
)
fit_4  = sampling(stan_model_4,
                  data=data,
                  chains=4,
                  iter=30000,
                  warmup=1000
)


posterior_1 = as.matrix(fit_1)
posterior_2 = as.matrix(fit_2)
posterior_3 = as.matrix(fit_3)
posterior_4 = as.matrix(fit_4)

density_1=density(posterior_1[, "b[2]"])
density_2=density(posterior_2[, "b[2]"])
density_3=density(posterior_3[, "b[2]"])
density_4=density(posterior_4[, "b[2]"])

png("posterior_change_prior.png", width=800, height=600, res=150)

par(mfrow=c(2,2), oma=rep(0.1, 4), mar=rep(1,4))

x=seq(-10,10, by=0.01)

plot(density_1, col="purple", xlim=c(-1,2), ylim=c(0,2.5), lwd=2, main="flat prior", xaxt="n", yaxt="n")
lines(density_1, col="red", xlim=c(-1,2), lwd=2)
lines(x, dnorm(x, mean=0, sd=10), col="blue", lwd=2)

legend("topleft",legend=c("prior", "likelihood", "posterior"), bty="n", lty=c(1,1,1), lwd=c(2,2,2), col=c("blue","red","purple"))

plot(density_2, col="purple", xlim=c(-1,2), ylim=c(0,2.5), lwd=2, main="weakly informative prior", xaxt="n", yaxt="n")
lines(density_1, col="red", xlim=c(-1,2), lwd=2)
lines(x, dnorm(x, mean=0, sd=1.2), col="blue", lwd=2)

plot(density_3, col="purple", xlim=c(-1,2), ylim=c(0,2.5), lwd=2, main="informative prior", xaxt="n", yaxt="n")
lines(density_1, col="red", xlim=c(-1,2), lwd=2)
lines(x, dnorm(x, mean=0, sd=0.5), col="blue", lwd=2)

plot(density_4, col="purple", xlim=c(-1,2), ylim=c(0,2.5), lwd=2, main="highly informative prior", xaxt="n", yaxt="n")
lines(density_1, col="red", xlim=c(-1,2), lwd=2)
lines(x, dnorm(x, mean=0, sd=0.2), col="blue", lwd=2)

dev.off()

#------------------------------------------------------------------------------

n=20
x = runif(n=n, min=-1, max=1)
y = rnorm(n=n, mean=a+b*x, sd=sigma)
df = data.frame(x=x,
                y=y)
data = list(n=n, 
            x=df$x, 
            y=df$y)
fit_1  = sampling(stan_model_3,
                  data=data,
                  chains=4,
                  iter=30000,
                  warmup=1000
)
fit_1_flat  = sampling(stan_model_1,
                  data=data,
                  chains=4,
                  iter=30000,
                  warmup=1000
)


n=50
x = runif(n=n, min=-1, max=1)
y = rnorm(n=n, mean=a+b*x, sd=sigma)
df = data.frame(x=x,
                y=y)
data = list(n=n, 
            x=df$x, 
            y=df$y)
fit_2  = sampling(stan_model_3,
                  data=data,
                  chains=4,
                  iter=30000,
                  warmup=1000
)
fit_2_flat  = sampling(stan_model_1,
                       data=data,
                       chains=4,
                       iter=30000,
                       warmup=1000
)


n=100
x = runif(n=n, min=-1, max=1)
y = rnorm(n=n, mean=a+b*x, sd=sigma)
df = data.frame(x=x,
                y=y)
data = list(n=n, 
            x=df$x, 
            y=df$y)
fit_3  = sampling(stan_model_3,
                  data=data,
                  chains=4,
                  iter=30000,
                  warmup=1000
)
fit_3_flat  = sampling(stan_model_1,
                       data=data,
                       chains=4,
                       iter=30000,
                       warmup=1000
)


n=500
x = runif(n=n, min=-1, max=1)
y = rnorm(n=n, mean=a+b*x, sd=sigma)
df = data.frame(x=x,
                y=y)
data = list(n=n, 
            x=df$x, 
            y=df$y)
fit_4  = sampling(stan_model_3,
                  data=data,
                  chains=4,
                  iter=30000,
                  warmup=1000
)
fit_4_flat  = sampling(stan_model_1,
                       data=data,
                       chains=4,
                       iter=30000,
                       warmup=1000
)


posterior_1 = as.matrix(fit_1)
posterior_2 = as.matrix(fit_2)
posterior_3 = as.matrix(fit_3)
posterior_4 = as.matrix(fit_4)

posterior_1_flat = as.matrix(fit_1_flat)
posterior_2_flat = as.matrix(fit_2_flat)
posterior_3_flat = as.matrix(fit_3_flat)
posterior_4_flat = as.matrix(fit_4_flat)

density_1=density(posterior_1[, "b[2]"])
density_2=density(posterior_2[, "b[2]"])
density_3=density(posterior_3[, "b[2]"])
density_4=density(posterior_4[, "b[2]"])

density_1_flat=density(posterior_1_flat[, "b[2]"])
density_2_flat=density(posterior_2_flat[, "b[2]"])
density_3_flat=density(posterior_3_flat[, "b[2]"])
density_4_flat=density(posterior_4_flat[, "b[2]"])

png("posterior_change_n.png", width=800, height=600, res=150)

par(mfrow=c(2,2), oma=rep(0.1, 4), mar=rep(1,4))

x=seq(-10,10, by=0.01)

plot(density_1, col="purple", xlim=c(-1,2), ylim=c(0,2.5), lwd=2, main="n=10", xaxt="n", yaxt="n")
lines(density_1_flat, col="red", lwd=2)
lines(x, dnorm(x, mean=0, sd=0.5), col="blue", lwd=2)

legend("topleft",legend=c("prior", "likelihood", "posterior"), bty="n", lty=c(1,1,1), lwd=c(2,2,2), col=c("blue","red","purple"))

plot(density_2, col="purple", xlim=c(-1,2), ylim=c(0,2.5), lwd=2, main="n=50", xaxt="n", yaxt="n")
lines(density_2_flat, col="red", lwd=2)
lines(x, dnorm(x, mean=0, sd=0.5), col="blue", lwd=2)

plot(density_3, col="purple", xlim=c(-1,2), ylim=c(0,2.5), lwd=2, main="n=100", xaxt="n", yaxt="n")
lines(density_3_flat, col="red", lwd=2)
lines(x, dnorm(x, mean=0, sd=0.5), col="blue", lwd=2)

plot(density_4, col="purple", xlim=c(-1,2), ylim=c(0,2.5), lwd=2, main="n=500", xaxt="n", yaxt="n")
lines(density_4_flat, col="red", lwd=2)
lines(x, dnorm(x, mean=0, sd=0.5), col="blue", lwd=2)

dev.off()