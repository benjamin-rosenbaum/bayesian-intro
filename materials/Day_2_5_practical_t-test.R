rm(list=ls())
library(rstan)
library(coda)
library(BayesianTools)
setwd("~/Desktop/teaching Bayes")
 
rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

# In this exercise, we will learn different ways to code a categorical predictor. 
# We will also learn how to use the `generated quantities{}` block.

#------------------------------------------------------------------------------
# statistical model 
#------------------------------------------------------------------------------
# Suppose we have measurements of a continuous variable in 2 groups,
# e.g. individual body masses in 2 populations.
# We want to test if both groups have a different mean.
# statistical model:
# y1_i ~ normal(mu_1, sigma) i=1,...,n1
# y2_i ~ normal(mu_2, sigma) i=1,...,n2
# notice that here, we assume that both populations have the same standard deviation.
# we could also model 2 separate standard deviations.

# Research question: what is the mean difference $\delta=\mu_2-\mu_1$ of both populations?

#------------------------------------------------------------------------------
# generate data 
#------------------------------------------------------------------------------
set.seed(123) # initiate random number generator for reproducability

n.1 = 30
mu.1 = 1
sigma.1 = 1

n.2 = 40
mu.2 = 2
sigma.2 = 1

y.1 = rnorm(n=n.1, mean=mu.1, sd=sigma.1)
y.2 = rnorm(n=n.2, mean=mu.2, sd=sigma.2)

par(mfrow=c(1,2))
hist(y.1)
hist(y.2)

#------------------------------------------------------------------------------
# straightforward model
#------------------------------------------------------------------------------

# We translate the statistical model above straightforward into a Stan model.
# Diffence in means will be investigated after fitting.

data = list(n1=n.1, 
            n2=n.2, 
            y1=y.1, 
            y2=y.2)

stan_code = '
data {
  int n1;
  vector[n1] y1;
  int n2;
  vector[n2] y2;
}
parameters {
  real mu1;  
  real mu2;  
  real<lower=0> sigma1;
  real<lower=0> sigma2;
}
model {
  // priors
  mu1 ~ normal(0, 10);
  mu2 ~ normal(0, 10);
  sigma1 ~ normal(0, 10);
  sigma2 ~ normal(0, 10);
  // likelihood
  y1 ~ normal(mu1, sigma1);
  y2 ~ normal(mu2, sigma2);
}
'

stan_model = stan_model(model_code=stan_code)
# save(stan_model, file="stan_model.Rdata")
# load("stan_model.Rdata")

fit  = sampling(stan_model,
                data=data,
                chains=3,
                iter=2000,
                warmup=1000
)

# Model convergence diagnostics look OK!

print(fit, digits=3, probs=c(0.025, 0.975))

plot(fit)
plot(As.mcmc.list(fit)) # from coda package


# To answer our research questions, we look at the posterior distribution 
# of $\delta=\mu_2-\mu_1$.
# First, we convert the posterior into a matrix. 
# Each row represents a sample from the posterior distribution.

posterior=as.matrix(fit)

str(posterior)
head(posterior)

correlationPlot(posterior[, 1:3], thin=1)

# Now we compute the distribution of differences. 
# The estimated mean difference in populations is 1.20 
# with a standard deviation (uncertainty) of 0.22.
# We compute the posterior probability of \delta being postitive given the data P(\delta>0|y): 
# We count the number the event occurs and divide by the number of samples. 
# Given the data, we are 100% sure that \delta>0, i.e. \mu_1 < \mu_2

delta = posterior[, "mu2"] - posterior[, "mu1"]

mean(delta)
sd(delta)

hist(delta)
abline(v=0, col="red", lwd=2)
sum(delta>0)/nrow(posterior)

#------------------------------------------------------------------------------
# same stan model, use "generated quantities" for calculating the difference
#------------------------------------------------------------------------------

# Instead of computing delta after fitting, we can include it in the model.

stan_code_delta = '
data {
  int n1;
  vector[n1] y1;
  int n2;
  vector[n2] y2;
}
parameters {
  real mu1;  
  real mu2;  
  real<lower=0> sigma;  
}
model {
  // priors
  mu1 ~ normal(0, 10);
  mu2 ~ normal(0, 10);
  sigma1 ~ normal(0, 10);
  sigma2 ~ normal(0, 10);
  // likelihood
  y1 ~ normal(mu1, sigma1);
  y2 ~ normal(mu2, sigma2);
}
generated quantities{
  real delta;
  delta = mu2-mu1;
}
'

stan_model_delta = stan_model(model_code=stan_code_delta)
# save(stan_model_delta, file="stan_model_delta.Rdata")
# load("stan_model_delta.Rdata")

fit_delta  = sampling(stan_model_delta,
                      data=data,
                      chains=3,
                      iter=2000,
                      warmup=1000
)

# Now we can directly look at delta:

print(fit_delta, digits=3, probs=c(0.025, 0.975))

plot(fit_delta)
plot(As.mcmc.list(fit_delta)) # from coda package

#------------------------------------------------------------------------------
# "classical" data format (as for lm)
#------------------------------------------------------------------------------

# Here we code the data differently, just one dataset with a factorial predictor "group"

df = data.frame(y=c(y.1,y.2), group=c(rep(1,n.1),rep(2,n.2)))

head(df)

data = list(y=df$y,
            group=df$group,
            n=nrow(df))

# Note that the factorial group is an integer (1 or 2), 
# so we can use it as an index in the Stan model.

# The statistical model reads:
# y_i ~ normal(mu_group_i, sigma)
# mu now is a vector of length 2.

stan_code_alternative = '
data {
  int n;
  real y[n]; // vector[n] y;
  int group[n];
}
parameters {
  real mu[2];  
  real<lower=0> sigma;  
}
model {
  // priors
  mu ~ normal(0, 10);
  sigma ~ normal(0, 10);
  // likelihood
  for(i in 1:n){
    y[i] ~ normal(mu[ group[i] ] , sigma);
  }
}
generated quantities{
  real delta;
  delta = mu[2]-mu[1];
}
'

stan_model_alternative = stan_model(model_code=stan_code_alternative)
# save(stan_model_alternative, file="stan_model_alternative.Rdata")
# load("stan_model_alternative.Rdata")

fit_alternative  = sampling(stan_model_alternative,
                            data=data,
                            chains=3,
                            iter=2000,
                            warmup=1000
)

print(fit_alternative, digits=3, probs=c(0.025, 0.975))

#------------------------------------------------------------------------------
# "classical" data format, "classical" model formulation with dummy coding
#------------------------------------------------------------------------------

# Again, we code the data a little differently. 
# Same as above, but group=0 for population 1 and group=1 for population 2.
# This is known as dummy coding and the statistical model reads
# y_i ~ normal(mu + delta*group_i, sigma)
# Here, mu is a single value and decribes the mean of population 1 (group=0),
# delta is the "effect" of population 2 compared to population 1, 
# the mean of population 2 is mu + delta (group=1)

df = data.frame(y=c(y.1,y.2), group=c(rep(0,n.1),rep(1,n.2))) # now group1=0, group2=1

head(df)

data = list(y=df$y,
            group=df$group,
            n=nrow(df))

stan_code_alternative_2 = '
data {
  int n;
  real y[n];
  int group[n];
}
parameters {
  real mu;  
  real delta;
  real<lower=0> sigma;  
}
model {
  // priors
  mu ~ normal(0, 10);
  delta ~ normal(0, 10);
  sigma ~ normal(0, 10);
  // likelihood
  for(i in 1:n){
    y[i] ~ normal(mu + group[i]*delta, sigma);  
  }
}
'

stan_model_alternative_2 = stan_model(model_code=stan_code_alternative_2)
# save(stan_model_alternative_2, file="stan_model_alternative_2.Rdata")
# load("stan_model_alternative_2.Rdata")

fit_alternative_2  = sampling(stan_model_alternative_2,
                              data=data,
                              chains=3,
                              iter=2000,
                              warmup=1000
)

print(fit_alternative_2, digits=3, probs=c(0.025, 0.975))

#------------------------------------------------------------------------------
# frequentist solutions
#------------------------------------------------------------------------------

t.test(y.1, y.2)

df$group = as.factor(df$group)
summary(lm(y~group, data=df))

# Common statistical tests are linear models!
# https://lindeloev.github.io/tests-as-linear/
