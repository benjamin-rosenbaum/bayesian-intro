rm(list=ls())
library(rstan)
library(coda)
library(BayesianTools)
library(boot)
setwd("~/Desktop/teaching Bayes")

set.seed(123) # initiate random number generator for reproducability

rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

#------------------------------------------------------------------------------
# statistical model 
#------------------------------------------------------------------------------
# Suppose we want to measure the effect of a continuous variable on presence / absence data.
# E.g. presence of a species in different locations y\in[0,1] presence (=1) or absence (=0), x = temperature in location.

# Statistical model: 
# y_i ~ bernoulli(p_i) = binomial(1,p_i)
# logit(p_i) = a+b*x_i
# or, equivalently: p_i = inverse.logit(a+b*x_i)

# p_i is the probability of presence of the species, the distribution of y_i is bernoulli (coin-flip).
# We assume a linear relationship of logit(p) with temperature. (Realistically, a hump-shaped relationship would make more sense.)
# inverse.logit transforms the values of the whole axis to the interval (0,1)

# Research question: is there a positive effect of the predictor to the presence of the species?

par(mfrow=c(1,2))
curve(logit, from=0, to=1)
curve(inv.logit, from=-10, to=10)

#------------------------------------------------------------------------------
# generate data 
#------------------------------------------------------------------------------
n = 50
x = sort(runif(n, 0, 1))

a = -2
b = 6

p = inv.logit(a+b*x)

y = rbinom(n=n, size=1, prob=p)

par(mfrow=c(1,1))
plot(x,y)
lines(x,p, lty=2, col="red")

#------------------------------------------------------------------------------
# straightforward model
#------------------------------------------------------------------------------
data = list(n=n, 
            x=x, 
            y=y)

stan_code = '
data {
  int  n;
  real x[n];
  int  y[n];
}
parameters {
  real a;  
  real b;  
}
model {
  real p[n];
  // priors
  a ~ normal(0, 10);
  b ~ normal(0, 10);
  // likelihood
  for(i in 1:n){
    p[i] = inv_logit( a+b*x[i] );
    y[i] ~ bernoulli(p[i]);
  }
}
'

stan_model = stan_model(model_code=stan_code)
# save(stan_model, file="stan_model.RData")
# load("stan_model.RData")

fit  = sampling(stan_model,
                data=data,
                chains=3,
                iter=3000,
                warmup=1000
)

print(fit, digits=3, probs=c(0.025, 0.975))

plot(fit)
plot(As.mcmc.list(fit)) # from coda package

posterior=as.matrix(fit)

head(posterior)

correlationPlot(posterior[, 1:2], thin=1)

#------------------------------------------------------------------------------
# predictions 
#------------------------------------------------------------------------------

# First, we generate credible intervals for the determinstic model (for p). 
# (90%, but choose as you like)
# later, we generate prediction intervals for the data (for y) 
# using also the stochastic part.

x.pred = seq(from=0, to=1, by=0.01)
p.pred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))

for(i in 1:nrow(posterior)){
  p.pred[i, ] = inv.logit(posterior[i,"a"] + posterior[i,"b"]*x.pred)
}

plot(x,y)
for(i in 1:100){
  lines(x.pred, p.pred[i, ], col=adjustcolor("red", alpha.f=0.3))
}

plot(x,y, xlim=c(0,1))

p.pred.mean = apply(p.pred, 2, function(x) mean(x)) 
lines(x.pred, p.pred.mean, col="red", lwd=2)

p.pred.q05 = apply(p.pred, 2, function(x) quantile(x, probs=0.05)) 
lines(x.pred, p.pred.q05, col="red", lwd=2, lty=2)

p.pred.q95 = apply(p.pred, 2, function(x) quantile(x, probs=0.95)) 
lines(x.pred, p.pred.q95, col="red", lwd=2, lty=2)

# lines(x.pred, inv.logit(a+b*x.pred), lwd=2, lty=1, col="blue")

# Now, we draw predicted data from the bernoulli (=binomial(1)) distribution in the statistical model.

y.pred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))
for(i in 1:nrow(posterior)){
  y.pred[i, ] = rbinom(n=length(x.pred), size=1, p=inv.logit(posterior[i,"a"] + posterior[i,"b"]*x.pred) )
  # y.pred[i, ] = rbinom(n=length(x.pred), size=1, p=p.pred[i, ])
}

y.pred.mean = apply(y.pred, 2, function(x) mean(x)) 
lines(x.pred, y.pred.mean, col="blue", lwd=2)

y.pred.q05 = apply(y.pred, 2, function(x) quantile(x, probs=0.05)) 
lines(x.pred, y.pred.q05, col="blue", lwd=2, lty=2)

y.pred.q95 = apply(y.pred, 2, function(x) quantile(x, probs=0.95)) 
lines(x.pred, y.pred.q95, col="blue", lwd=2, lty=2)

#------------------------------------------------------------------------------
# vectorised stan code
#------------------------------------------------------------------------------
stan_code_2 = '
data {
  int n;
  vector[n] x;
  int y[n];
}
parameters {
  real a;  
  real b;  
}
model {
  // priors
  a ~ normal(0, 10);
  b ~ normal(0, 10);
  // likelihood
  y ~ bernoulli(inv_logit(a+b*x));
}
'

stan_model_2 = stan_model(model_code=stan_code_2)
# save(stan_model_2, file="stan_model_2.RData")
# load("stan_model_2.RData")

fit_2  = sampling(stan_model_2,
                  data=data,
                  chains=3,
                  iter=2000,
                  warmup=1000
)

print(fit_2)


#------------------------------------------------------------------------------
# frequentist solution
#------------------------------------------------------------------------------
summary(glm(y~x, family=binomial(link="logit")))

