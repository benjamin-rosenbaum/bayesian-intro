rm(list=ls())
library(rstan)
library(coda)
library(BayesianTools)
library(boot)
setwd("~/Desktop/teaching Bayes")
 
rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

# In this exercise, we move to a nonlinear model in population dynamics. 
# The deterministic model may look complicated, but the statistics are the same as before.

curve(1*exp(2*x), xlab="t", ylab="y(t)", main="exponential growth")

# Exponential growth of a population density y(t) over time t is described by a differential equation
# dy(t)/dt = r*y(t),
# with an initial density y(t0)=y0 in t0=0 and a growth rate r. 
# There is an explicit formulation that describes the population density at time t:
# y(t)=y0*exp(r*t)
# This model is unrealistic in the sense that it grows exponentially to infinity over time.

curve(10/( 1+(10-1)/1 * exp(-2*x) ) , from=0, to=4, xlab="t", ylab="y(t)", main="logistic growth")

# A more realistic model is logistic growth, where the population densitiy grows only it reaches a carrying capacity K.
# dy(t)/dt = r*(1-y(t)/K)*y(t),
# with an initial density y(t0)=y0 in t0=0. If carrying capacity is reached (y(t)=K), 
# the term in parentheses and therefore the whole expression is zero and the population size reaches an equilibrium.
# For this differential equation, there also is an explicit formulation describing population size y(t) at time t:
# y(t) = K/( 1+(K-y0)/y0 * exp(-r*t) )
# This looks complicated, but for any given y0,K,r we can directly compute y(t).

# Given some time series data t_i, y_i (i=1,...,n), we formulate a statistical model:
# mu_i = K/( 1+(K-y0)/y0 * exp(-r*t_i) )
# y_i ~ normal(mu_i, sigma_i)
# The parameters are y0, r, K. 
# Why is y0 a parameter? There is measurement for that (observation in t0). 
# Yes, ans that measurement is assumed to have an observation error sigma just like all observations!
# Have a look at the equation, y0 has an influence on all y(t), so leaving it as a free parameter leaves more freedom to the model.
# Just think of it as an intercept in t_0!

# There are more sophisticated ways to model time series data, but for now let's just fit this nonlinear equation y(t) with a predictor t.


#------------------------------------------------------------------------------
# generate data 
#------------------------------------------------------------------------------
set.seed(123) # initiate random number generator for reproducability

n = 50
t = 0:(n-1)

r = 0.2
K = 6
y0 = 0.1

y = K/( 1+(K-y0)/y0 * exp(-r*t) )

y = exp(rnorm(n,mean=log(y),sd=0.1))

plot(t,y)

data = list(n=n,
            y=y,
            t=t)

#------------------------------------------------------------------------------
# straightforward model
#------------------------------------------------------------------------------
stan_code = '
data {
  int n;
  real y[n];
  real t[n];
}
parameters {
  real<lower=0> y0;  
  real<lower=0> r;
  real<lower=0> K;
  real<lower=0> sigma;
}
model {
  // priors
  y0 ~ normal(0, 10);
  r ~ normal(0, 10);
  K ~ normal(0, 10);
  sigma ~ normal(0,10);

  // likelihood
  for (i in 1:n){
    y[i] ~ normal( K/( 1+(K-y0)/y0 * exp(-r*t[i]) ), sigma );
  }
}
'

stan_model = stan_model(model_code=stan_code)
# save(stan_model, file="stan_model.RData")
# load("stan_model.RData")

fit  = sampling(stan_model,
                data=data,
                chains=3,
                iter=2000,
                warmup=1000
)

print(fit, digits=3, probs=c(0.025, 0.975))

# plot(fit)
plot(As.mcmc.list(fit)) # from coda package

posterior=as.matrix(fit)

correlationPlot(posterior[,1:4], thin=1)


#------------------------------------------------------------------------------
# predictions 
#------------------------------------------------------------------------------

# We generate credible / confidence intervals for the deterministic model and prediction intervals for the whole model.
# Again, you could also code the directly in the model using the "generated quantitites{}" block.

t.pred = data$t
y.cred = matrix(0, nrow=nrow(posterior), ncol=length(t.pred))

for(i in 1:nrow(posterior)){
  y.cred[i, ] = posterior[i,"K"]/( 1+(posterior[i,"K"]-posterior[i,"y0"])/posterior[i,"y0"] * exp(-posterior[i,"r"]*t.pred) )
}

plot(t,y)

y.cred.mean = apply(y.cred, 2, function(x) mean(x)) 
lines(t.pred, y.cred.mean, col="red", lwd=2)

y.cred.q05 = apply(y.cred, 2, function(x) quantile(x, probs=0.05)) 
lines(t.pred, y.cred.q05, col="red", lwd=2, lty=2)

y.cred.q95 = apply(y.cred, 2, function(x) quantile(x, probs=0.95)) 
lines(t.pred, y.cred.q95, col="red", lwd=2, lty=2)

y.pred = matrix(0, nrow=nrow(posterior), ncol=length(t.pred))
for(i in 1:nrow(posterior)){
  y.pred[i, ] = rnorm(n=length(t.pred), 
                     mean = posterior[i,"K"]/( 1+(posterior[i,"K"]-posterior[i,"y0"])/posterior[i,"y0"] * exp(-posterior[i,"r"]*t.pred) ),
                     sd = posterior[i,"sigma"])
}

y.pred.mean = apply(y.pred, 2, function(x) mean(x)) 
lines(t.pred, y.pred.mean, col="blue", lwd=2)

y.pred.q05 = apply(y.pred, 2, function(x) quantile(x, probs=0.05)) 
lines(t.pred, y.pred.q05, col="blue", lwd=2, lty=2)

y.pred.q95 = apply(y.pred, 2, function(x) quantile(x, probs=0.95)) 
lines(t.pred, y.pred.q95, col="blue", lwd=2, lty=2)

# The model predicts values < 0. For populations, this doesn't make any sense. 

# Let's look closer at observed vs. predicted

# here, t.pred are the same predictor values as in the data. so we can make the observed-predicted plot directly.

plot(data$y, y.pred.mean, ylim=c(min(data$y), max(data$y)))
abline(0,1)
for (i in 1:n){
  lines(c(data$y[i], data$y[i]), c(y.pred.q05[i], y.pred.q95[i]))
}

# stochastic model: y_i~normal(mu_i,sigma)
# Even if the deterministic model ensures mu_i>0, the statistical model (normal distribution) allows negative predictions, while it should be bounded at zero.
# Also, the variation in the data seems to grow with y ("heteroscedasticity").

# Next, we will try a more appropriate stochastic model (lognormal distribution)

# y_i ~ lognormal( log(\mu_i), sigma ), or equivalently
# log(y_i) ~ normal( log(\mu_i), sigma ),
# so the logs of the (nonnegative) observations are normally distibuted. 
# Note that mean and standard deviation are now defined on the log-scale

#------------------------------------------------------------------------------
# lognormal residual model
#------------------------------------------------------------------------------
stan_code_lognormal = '
data {
  int n;
  real y[n];
  real t[n];
}
parameters {
  real<lower=0> y0;  
  real<lower=0> r;
  real<lower=0> K;
  real<lower=0> sigma;
}
model {
  // priors
  y0 ~ normal(0, 10);
  r ~ normal(0, 10);
  K ~ normal(0, 10);
  sigma ~ normal(0,10);
  
  // likelihood
  for (i in 1:n){
    y[i] ~ lognormal( log( K/( 1+(K-y0)/y0 * exp(-r*t[i]) ) ), sigma );
  }
}
'

stan_model_lognormal = stan_model(model_code=stan_code_lognormal)
# save(stan_model_lognormal, file="stan_model_lognormal.RData")
# load("stan_model_lognormal.RData")

fit_lognormal  = sampling(stan_model_lognormal,
                          data=data,
                          chains=3,
                          iter=2000,
                          warmup=1000
)

print(fit_lognormal, digits=3, probs=c(0.025, 0.975))

# plot(fit_lognormal)
plot(As.mcmc.list(fit_lognormal)) # from coda package

posterior=as.matrix(fit_lognormal)

#------------------------------------------------------------------------------
# predictions 
#------------------------------------------------------------------------------

# Again, credible / confidence intervals for the deterministic model and prediction intervals for the whole model.

t.pred = data$t
y.cred = matrix(0, nrow=nrow(posterior), ncol=length(t.pred))

for(i in 1:nrow(posterior)){
  y.cred[i, ] = posterior[i,"K"]/( 1+(posterior[i,"K"]-posterior[i,"y0"])/posterior[i,"y0"] * exp(-posterior[i,"r"]*t.pred) )
}

plot(t,y)

y.cred.mean = apply(y.cred, 2, function(x) mean(x)) 
lines(t.pred, y.cred.mean, col="red", lwd=2)

y.cred.q05 = apply(y.cred, 2, function(x) quantile(x, probs=0.05)) 
lines(t.pred, y.cred.q05, col="red", lwd=2, lty=2)

y.cred.q95 = apply(y.cred, 2, function(x) quantile(x, probs=0.95)) 
lines(t.pred, y.cred.q95, col="red", lwd=2, lty=2)

y.pred = matrix(0, nrow=nrow(posterior), ncol=length(t.pred))
for(i in 1:nrow(posterior)){
  y.pred[i, ] = rlnorm(n=length(t.pred), 
                      meanlog = log( posterior[i,"K"]/( 1+(posterior[i,"K"]-posterior[i,"y0"])/posterior[i,"y0"] * exp(-posterior[i,"r"]*t.pred) ) ),
                      sdlog = posterior[i,"sigma"])
}

y.pred.mean = apply(y.pred, 2, function(x) mean(x)) 
lines(t.pred, y.pred.mean, col="blue", lwd=2)

y.pred.q05 = apply(y.pred, 2, function(x) quantile(x, probs=0.05)) 
lines(t.pred, y.pred.q05, col="blue", lwd=2, lty=2)

y.pred.q95 = apply(y.pred, 2, function(x) quantile(x, probs=0.95)) 
lines(t.pred, y.pred.q95, col="blue", lwd=2, lty=2)

plot(data$y, y.pred.mean, ylim=c(min(data$y), max(data$y)))
abline(0,1)
for (i in 1:n){
  lines(c(data$y[i], data$y[i]), c(y.pred.q05[i], y.pred.q95[i]))
}

# Model predictions look much better than before. The variation in the predicted values grows with y, just as in the observed data.
# On the logscale, variation is constant:

plot(t,y, log="y")
lines(t.pred, y.cred.mean, col="red", lwd=2)
lines(t.pred, y.cred.q05, col="red", lwd=2, lty=2)
lines(t.pred, y.cred.q95, col="red", lwd=2, lty=2)

lines(t.pred, y.pred.mean, col="blue", lwd=2)
lines(t.pred, y.pred.q05, col="blue", lwd=2, lty=2)
lines(t.pred, y.pred.q95, col="blue", lwd=2, lty=2)

plot(data$y, y.pred.mean, ylim=c(min(data$y), max(data$y)), log="xy")
abline(0,1)
for (i in 1:n){
  lines(c(data$y[i], data$y[i]), c(y.pred.q05[i], y.pred.q95[i]))
}

