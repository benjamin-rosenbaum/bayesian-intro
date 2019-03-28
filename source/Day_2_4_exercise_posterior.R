rm(list=ls())
library(rstan)
library(coda)
library(BayesianTools)
setwd("~/Desktop/teaching Bayes")

# set options to use parallelization to run several chains simultaneously
rstan_options(auto_write = TRUE)
options(mc.cores = 3) # number of CPU cores

#------------------------------------------------------------------------------
# generate data 
#------------------------------------------------------------------------------
set.seed(123) # initiate random number generator for reproducability

n=50

a=1.0
b=0.5
c=0.4
sigma=0.15

x = runif(n=n, min=0, max=1)
y = rnorm(n=n, mean=a+b*x+c*x^2, sd=sigma)

df = data.frame(x=x,
                y=y)

plot(df)

#------------------------------------------------------------------------------
# 2 statistical model, deterministic and stochastic part
#------------------------------------------------------------------------------

# statistical model: linear regression
# y_i ~ normal(\mu_i, \sigma)
# mu_i = a+b*x_i+c*x_i^2

stan_code_lin = '
data {
  int n;
  vector[n] x;
  vector[n] y;
}
parameters {
  real a;  
  real b;  
  real<lower=0> sigma;  // standard deviation
}
model {
  // priors
  a ~ normal(0, 10);
  b ~ normal(0, 10);
  sigma ~ normal(0, 10);
  // likelihood
  y ~ normal(a+b*x, sigma);
}
'

stan_code_quad = '
data {
  int n;
  vector[n] x;
  vector[n] y;
}
parameters {
  real a;  
  real b;  
  real c;
  real<lower=0> sigma;  // standard deviation
}
model {
  // priors
  a ~ normal(0, 10);
  b ~ normal(0, 10);
  c ~ normal(0, 10);
  sigma ~ normal(0, 10);
  // likelihood
  y ~ normal(a + b*x + c * x .* x , sigma);
}
'

#------------------------------------------------------------------------------
# prepare data and stan options
#------------------------------------------------------------------------------
data = list(n=n, 
            x=df$x, 
            y=df$y)


stan_model_lin = stan_model(model_code=stan_code_lin)
stan_model_quad = stan_model(model_code=stan_code_quad)

#------------------------------------------------------------------------------
# MCMC sampling 
#------------------------------------------------------------------------------
# fit.1  = sampling(stan_model_lin,
#                 data=data,
#                 chains=3,
#                 iter=2000,
#                 warmup=1000
# )

fit.2  = sampling(stan_model_quad,
                data=data,
                chains=3,
                iter=2000,
                warmup=1000
)

#------------------------------------------------------------------------------
# explore the posterior 
#------------------------------------------------------------------------------
print(fit.2, digits=3, probs=c(0.025, 0.975))

plot(fit.2)
plot(As.mcmc.list(fit.2)) # from coda package

posterior=as.matrix(fit.2)

correlationPlot(posterior[, 1:4], thin=1) # from BayesianTools package

hist(posterior[, "c"])
abline(v=0, col="red", lwd=2)
sum(posterior[, "c"]>0)/nrow(posterior)

#------------------------------------------------------------------------------
# posterior predictions 
#------------------------------------------------------------------------------
x.pred = seq(from=0, to=1, by=0.1)
y.cred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))

for(i in 1:nrow(posterior)){
  y.cred[i, ] = posterior[i,"a"] + posterior[i,"b"]*x.pred + posterior[i,"c"]*x.pred^2
}

plot(df)
for(i in 1:100){
  lines(x.pred, y.cred[i, ], col=adjustcolor("red", alpha.f=0.3))
}

plot(df)

y.cred.mean = apply(y.cred, 2, function(x) mean(x)) 
lines(x.pred, y.cred.mean, col="red", lwd=2)

y.cred.q05 = apply(y.cred, 2, function(x) quantile(x, probs=0.05)) 
lines(x.pred, y.cred.q05, col="red", lwd=2, lty=2)

y.cred.q95 = apply(y.cred, 2, function(x) quantile(x, probs=0.95)) 
lines(x.pred, y.cred.q95, col="red", lwd=2, lty=2)

y.pred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))
for(i in 1:nrow(posterior)){
  y.pred[i, ] = rnorm(n=length(x.pred), mean=y.cred[i, ], sd=rep(posterior[i, "sigma"],length(x.pred)) )
}

y.pred.mean = apply(y.pred, 2, function(x) mean(x)) 
lines(x.pred, y.pred.mean, col="blue", lwd=2)

y.pred.q05 = apply(y.pred, 2, function(x) quantile(x, probs=0.05)) 
lines(x.pred, y.pred.q05, col="blue", lwd=2, lty=2)

y.pred.q95 = apply(y.pred, 2, function(x) quantile(x, probs=0.95)) 
lines(x.pred, y.pred.q95, col="blue", lwd=2, lty=2)

#------------------------------------------------------------------------------
# posterior predictions - obs vs pred for data
#------------------------------------------------------------------------------

x.pred = df$x
y.cred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))

for(i in 1:nrow(posterior)){
  y.cred[i, ] = posterior[i,"a"] + posterior[i,"b"]*x.pred + posterior[i,"c"]*x.pred^2
}

y.cred.mean = apply(y.cred, 2, function(x) mean(x)) 
y.cred.q05 = apply(y.cred, 2, function(x) quantile(x, probs=0.05)) 
y.cred.q95 = apply(y.cred, 2, function(x) quantile(x, probs=0.95)) 

y.pred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))
for(i in 1:nrow(posterior)){
  y.pred[i, ] = rnorm(n=length(x.pred), mean=y.cred[i, ], sd=rep(posterior[i, "sigma"],length(x.pred)) )
}

y.pred.mean = apply(y.pred, 2, function(x) mean(x)) 
y.pred.q05 = apply(y.pred, 2, function(x) quantile(x, probs=0.05)) 
y.pred.q95 = apply(y.pred, 2, function(x) quantile(x, probs=0.95)) 

plot(df$y, y.pred.mean, ylim=c(min(df$y), max(df$y)), xlab="observed", ylab="predicted")
abline(0,1)
for (i in 1:n){
  lines(c(df$y[i], df$y[i]), c(y.pred.q05[i], y.pred.q95[i]))
}

