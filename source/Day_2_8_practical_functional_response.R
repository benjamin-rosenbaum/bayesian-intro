rm(list=ls())
library(rstan)
library(coda)
library(BayesianTools)
setwd("~/Desktop/teaching Bayes")

rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

# In feeding interactions between a predator and a prey species, the functional response F describes the 
# density-dependent feeding rate of the consumer. The number of prey individuals a predator can consume 
# in a given time does not grow proportional with increasing prey availability, but reaches a saturation.
# (A predator can only consume a limited number of prey in a given time.)
# Mathematically, F(N) is parameterised by attack rate a and handling time h
# F(N) = a*N/(1+a*h*N)
# a descibes the slope of the curve in the origin, while the saturation (maximum feeding rate) is given by 1/h.
# This is a nonlinear function and we can't use "lm()" to fit it.

a=1
h=0.1
curve(a*x/(1+a*h*x), from=0, to=100,
      xlab="N_0", ylab="N_eaten", main="type II functional response")

#------------------------------------------------------------------------------
# read data 
#------------------------------------------------------------------------------

# We read data from a feeding experiment. 
# Rosenbaum B, Rall BC. Fitting functional responses: Direct parameter estimation by simulating differential equations. Methods Ecol Evol. 2018;9:2076–2090. https://doi.org/10.1111/2041-210X.13039
# oringinal data from: Vucic‐Pestic, O. , Rall, B. C., Kalinkat, G. and Brose, U. (2010), Allometric functional response model: body masses constrain interaction strengths. Journal of Animal Ecology, 79: 249-256. doi:10.1111/j.1365-2656.2009.01622.x
# In repeated feeding trials with varying number of initial prey, the number of eaten prey was counted

df = read.csv("https://datadryad.org/bitstream/handle/10255/dryad.181841/data_vucic-pestic_manual.csv")
head(df)

plot(jitter(df[,1]), jitter(df[,2]),xlab="N_0", ylab="N_eaten")

# Let's assume that the experiment was supervised and each eaten prey item was replaced immediately, 
# such that the number of available prey was constant during the course of the experiment. 
# (Otherwise we would have to correct for prey depletion, see Rosenbaum & Rall, 2018).

# The deterministic model is defined by the functional response F(N). 
# This is count data, so we use a poisson distribution for the stochastic part.
# N_eaten_i ~ poisson(a*N_i/(1+a*h*N_i))

data = list(n = nrow(df),
            x = df$N0,
            y = df$Neaten)

stan_code = '
data {
  int n;
  int x[n];
  int y[n];
}
parameters {
  real<lower=0> a;  
  real<lower=0> h;  
}
model {
  // priors
  a ~ normal(0, 10);
  h ~ normal(0, 1);
  // likelihood
  for(i in 1:n){
    y[i] ~ poisson( a*x[i]/(1.0+a*h*x[i]) );
  }
}
'

stan_model = stan_model(model_code=stan_code)
save(stan_model, file="stan_model.RData")
# load("stan_model.RData")

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

head(posterior)

correlationPlot(posterior[,1:2], thin=1)


#------------------------------------------------------------------------------
# predictions 
#------------------------------------------------------------------------------

# 90% credible intervals for the deterministic part of the model

x.pred = seq(from=1, to=max(data$x), by=1)
y.cred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))

for(i in 1:nrow(posterior)){
  y.cred[i, ] = posterior[i,"a"]*x.pred / (1+posterior[i,"a"]*posterior[i,"h"]*x.pred)
}

plot(jitter(df[,1]), jitter(df[,2]))

for(i in 1:100){
  lines(x.pred, y.cred[i, ], col=adjustcolor("red", alpha.f=0.3))
}

plot(jitter(df[,1]), jitter(df[,2]))

y.cred.mean = apply(y.cred, 2, function(x) mean(x)) 
lines(x.pred, y.cred.mean, col="red", lwd=2)

y.cred.q05 = apply(y.cred, 2, function(x) quantile(x, probs=0.05)) 
lines(x.pred, y.cred.q05, col="red", lwd=2, lty=2)

y.cred.q95 = apply(y.cred, 2, function(x) quantile(x, probs=0.95)) 
lines(x.pred, y.cred.q95, col="red", lwd=2, lty=2)

# 90% predictions intervals for the data including stochastic part of the model

y.pred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))
for(i in 1:nrow(posterior)){
  y.pred[i, ] = rpois(n=length(x.pred), lambda=y.cred[i, ])
}

y.pred.mean = apply(y.pred, 2, function(x) mean(x)) 
lines(x.pred, y.pred.mean, col="blue", lwd=2)

y.pred.q05 = apply(y.pred, 2, function(x) quantile(x, probs=0.05)) 
lines(x.pred, y.pred.q05, col="blue", lwd=2, lty=2)

y.pred.q95 = apply(y.pred, 2, function(x) quantile(x, probs=0.95)) 
lines(x.pred, y.pred.q95, col="blue", lwd=2, lty=2)

# Variation in observed data increases with N, this variation is also captured in the predicted values.

