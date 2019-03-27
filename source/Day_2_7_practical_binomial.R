rm(list=ls())
library(rstan)
library(coda)
library(BayesianTools)
library(boot)
setwd("~/Desktop/teaching Bayes")


rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

#------------------------------------------------------------------------------
# statistical model 
#------------------------------------------------------------------------------
# Suppose we want to measure the effect of a continuous variable on presence / absence data.
# E.g. presence of a species in different locations, x = temperature in location.
# Now we have k trials per observation and y counts the number of successful trials, so y\in{0,1,2,...,k}

# Statistical model: 
# y_i ~ binomial(k,p_i)
# logit(p_i) = a+b*x_i
# or, equivalently: p_i = inverse.logit(a+b*x_i)

# p_i is the probability of presence of the species, the distribution of y_i is binomial with k trials.
# We assume a linear relationship of logit(p) with temperature. (Realistically, a hump-shaped relationship would make more sense.)
# inverse.logit transforms the values of the whole axis to the interval (0,1)

# Research question: is there a positive effect of the predictor to the presence of the species?

#------------------------------------------------------------------------------
# generate data 
#------------------------------------------------------------------------------
set.seed(123) # initiate random number generator for reproducability

n = 50
x = sort(runif(n, 0, 1))

a = -2
b = 5

p = inv.logit(a+b*x)

y = rbinom(n=n, size=4, prob=p)

plot(x,y)
lines(x,4*p)

#------------------------------------------------------------------------------
# straightforward model
#------------------------------------------------------------------------------
data = list(n=n, 
            x=x, 
            y=y)

stan_code = '
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
  y ~ binomial(4,inv_logit(a+b*x));
}
'

# stan_model = stan_model(model_code=stan_code)
# save(stan_model, file="stan_model.RData")
load("stan_model.RData")

fit  = sampling(stan_model,
                data=data,
                chains=3,
                iter=2000,
                warmup=1000
)

print(fit, digits=3, probs=c(0.025, 0.975))

plot(fit)
plot(As.mcmc.list(fit)) # from coda package

posterior=as.matrix(fit)

str(posterior)

#------------------------------------------------------------------------------
# predictions 
#------------------------------------------------------------------------------

# First, we generate credible intervals for the determinstic model. (90%, but choose as you like)
# The deterministic model is p = inv.logit(a+b*x), the probability of a successful trial.
# But there are k=4 trials, so the expected value of successful trials is 4*p.
# Later, we generate prediction intervals for the data (for y) using also the stochastic part.

x.pred = seq(from=0, to=1, by=0.01)
y.cred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))

for(i in 1:nrow(posterior)){
  y.cred[i, ] = 4 * inv.logit(posterior[i,"a"] + posterior[i,"b"]*x.pred)
}

plot(x,y)
for(i in 1:100){
  lines(x.pred, y.cred[i, ], col=adjustcolor("red", alpha.f=0.3))
}

plot(x,y)

y.cred.mean = apply(y.cred, 2, function(x) mean(x)) 
lines(x.pred, y.cred.mean, col="red", lwd=2)

y.cred.q05 = apply(y.cred, 2, function(x) quantile(x, probs=0.05)) 
lines(x.pred, y.cred.q05, col="red", lwd=2, lty=2)

y.cred.q95 = apply(y.cred, 2, function(x) quantile(x, probs=0.95)) 
lines(x.pred, y.cred.q95, col="red", lwd=2, lty=2)

# lines(x.pred, 4*inv.logit(a+b*x.pred), lwd=2, lty=1, col="blue")

# Now, we draw predicted data from the binomial distribution (k=4 trials) in the statistical model.

y.pred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))
for(i in 1:nrow(posterior)){
  y.pred[i, ] = rbinom(n=length(x.pred), size=4, p=inv.logit(posterior[i,"a"] + posterior[i,"b"]*x.pred))
}

y.pred.mean = apply(y.pred, 2, function(x) mean(x)) 
lines(x.pred, y.pred.mean, col="blue", lwd=2)

y.pred.q05 = apply(y.pred, 2, function(x) quantile(x, probs=0.05)) 
lines(x.pred, y.pred.q05, col="blue", lwd=2, lty=2)

y.pred.q95 = apply(y.pred, 2, function(x) quantile(x, probs=0.95)) 
lines(x.pred, y.pred.q95, col="blue", lwd=2, lty=2)

#------------------------------------------------------------------------------
# frequentist solution
#------------------------------------------------------------------------------

y.succ.fail = cbind(y, 4-y) # 1st column successes, 2nd column fails (if there are y successes, there must be k-y fails)
summary(glm(y.succ.fail~x, family=binomial(link="logit")))
