rm(list=ls())
library(rstan)
library(coda)
library(BayesianTools)
# setwd("~/Desktop/teaching Bayes")
 
set.seed(123) # initiate random number generator for reproducibility

rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

# In the previous model, we just fitted mean values to groups 
# and there was still some unexplained variation.
# Here, we will add a continuous predictor (covariate).
# Specifically, for the same dataset of fruitfly longevity and sexual activity, 
# we add the covariate individual body size.
# Lifespan is assumed to be positively correlated with body size.

#------------------------------------------------------------------------------
# load data 
#------------------------------------------------------------------------------

df = read.csv("~/git/bayesian-intro/data/FruitflyDataReduced.csv")
head(df)
str(df)

plot(df$Thorax, df$Longevity, col=as.factor(df$CompanionNumber))

# Note that the range of x-values (body size) is far away from zero. 
# This makes intercepts on regression hard to interpret.
# (what is the longevity of a fruitfly with zero body size?)
# We will use normalized (zscore) predictor values instead.

df$Thorax.norm = as.numeric(scale(df$Thorax))
df$group = as.integer(df$CompanionNumber)
plot(df$Thorax.norm, df$Longevity, col=as.factor(df$CompanionNumber))

# as.integer(<factor>) codes groups in alphabetical order
table(df$group, df$CompanionNumber)

par(mfrow=c(2,3))
for (i in 1:5){
  df.sub=subset(df, df$group==i)
  plot(df.sub$Thorax.norm,  df.sub$Longevity,
       xlim=range(df$Thorax.norm),
       ylim=range(df$Longevity),
       main=levels(df$CompanionNumber)[i]
  )
}

data = list(y = df$Longevity,
            x = df$Thorax.norm,
            group = df$group,
            n = nrow(df),
            n_group = 5)

#------------------------------------------------------------------------------
# partial pooling model 
#------------------------------------------------------------------------------

# We will fit linear regression lines to each group ("CompanionNumber) as follows:

# y_i ~ normal(a[group_i]+b*x_i, sigma) i=1,...,n (n observations)
# a_j ~ normal(mu_a, sigma_a) j=1,...,m (m groups)

# Here, a_j is a group-level intercept, which are allowed to vary (partial pooling). 
# (Analogue to mu_j in the previous example)
# But we assume identical slope b for all groups (complete pooling).
# So this is a random intercepts linear regression, 
# lm-formulation would be "y ~ (1|group) + x".
# The approach is also similar to frequentist ANCOVA.

# The Stan code differs from the previous model by adding 
# the covariate x and parameter slope b
# (and renaming mu to a, that's it!)

stan_code_partial = '
data {
  int n;
  int n_group;
  real y[n];
  real x[n];
  int group[n];
}
parameters {
  real a[n_group];  
  real b;  
  real<lower=0> sigma;
  real mu_a;
  real<lower=0> sigma_a;
}
model {
  // priors
  mu_a ~ normal(0,100);
  sigma_a ~ cauchy(0,10);
  for (j in 1:n_group){
    a[j] ~ normal(mu_a,sigma_a);
  }
  b ~ normal(0,100);
  sigma ~ normal(0,100);
  // likelihood
  for(i in 1:n){
    y[i] ~ normal(a[group[i]]+b*x[i], sigma); 
  }
}
'

stan_model_partial = stan_model(model_code=stan_code_partial)
# save(stan_model_partial, file="stan_code_partial.RData")
# load("stan_code_partial.RData")

fit_partial  = sampling(stan_model_partial,
                        data=data,
                        chains=3,
                        iter=2000,
                        warmup=1000
)

print(fit_partial, digits=3, probs=c(0.025, 0.975))

plot(fit_partial)
# plot(fit_partial, pars="a")
# plot(As.mcmc.list(fit_partial)) # from coda package

posterior=as.matrix(fit_partial)

# As before, we can look at the individual differences 
# of intercepts between groups ("contrasts").
# E.g., the posterior distribution of a4-a5

contrast = posterior[,"a[4]"]-posterior[,"a[5]"]

par(mfrow=c(1,1))

plot(density(contrast))
abline(v=0, col="red", lwd=2)

#------------------------------------------------------------------------------
# predictions / credible intervals
#------------------------------------------------------------------------------

# Again, we can generate predictions and compute credible intervals (for the deterministic part of the model)
# and prediction intervals (for the data, including statistical part of the model).
# (Here: 90% credible intervals)

x.pred = seq(from=min(df$Thorax.norm), to=max(df$Thorax.norm), by=0.1)

par(mfrow=c(2,3))

for (i in 1:5){
  df.sub=subset(df, df$group==i)
  plot(df.sub$Thorax.norm,  df.sub$Longevity,
       xlim=range(df$Thorax.norm),
       ylim=range(df$Longevity),
       main=levels(df$CompanionNumber)[i]
  )

  y.cred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))
  for(j in 1:nrow(posterior)){
    # column i in posterior corresponds to a_i, alternatively reference by name: 
    # posterior[j,paste0("a[",i,"]")]
    y.cred[j, ] = posterior[j,i] + posterior[j,"b"]*x.pred 
  }

  y.cred.mean = apply(y.cred, 2, function(x) mean(x)) 
  lines(x.pred, y.cred.mean, col="red", lwd=2)
  
  y.cred.q05 = apply(y.cred, 2, function(x) quantile(x, probs=0.05)) 
  lines(x.pred, y.cred.q05, col="red", lwd=2, lty=2)
  
  y.cred.q95 = apply(y.cred, 2, function(x) quantile(x, probs=0.95)) 
  lines(x.pred, y.cred.q95, col="red", lwd=2, lty=2)
}

