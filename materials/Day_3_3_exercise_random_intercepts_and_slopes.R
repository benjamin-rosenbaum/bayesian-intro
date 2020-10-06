rm(list=ls())
library(rstan)
library(coda)
library(BayesianTools)
# setwd("~/Desktop/teaching Bayes")

# remotes::install_github("allisonhorst/palmerpenguins")
library(palmerpenguins)

set.seed(123) # initiate random number generator for reproducability

rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

# Extend the last model by including random slopes, too!

#------------------------------------------------------------------------------
# load data 
#------------------------------------------------------------------------------

df = penguins
head(df)
str(df)

df = df[complete.cases(df), ]

df$body_mass_g_norm = as.numeric(scale(df$body_mass_g))
df$group = as.integer(df$species)

data = list(y = df$bill_length_mm,
            x = df$body_mass_g_norm,
            group = df$group,
            n = nrow(df),
            n_group = 3)

par(mfrow=c(3,1))
for (i in 1:3){
  df.sub=subset(df, df$group==i)
  plot(df.sub$body_mass_g_norm,  df.sub$bill_length_mm,
       xlim=range(df$body_mass_g_norm),
       ylim=range(df$bill_length_mm),
       main=levels(df$species)[i]
  )
}

#------------------------------------------------------------------------------
# partial pooling model 
#------------------------------------------------------------------------------

# We will fit linear regression lines to each group ("CompanionNumber) as follows:
# y_i ~ normal(a[group_i]+b[group_i]*x_i, sigma) i=1,...,n (n observations)
# a_j ~ normal(mu_a, sigma_a) j=1,...,m (m groups)
# b_j ~ normal(mu_b, sigma_b) j=1,...,m (m groups)
# Here, a_j and b_j are group-level intercepts and slopes,
# which are allowed to vary (partial pooling).
# Both have their own means and standard deviations, which are also free parameters to be estimated.
# So this is a random intercepts and slopes linear regression, lm-formulation would be `y ~ x + (x|group)`, 
# which is short for `y ~ 1+x + (1+x|group)`.

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
  real b[n_group];  
  real<lower=0> sigma;
  real mu_a;
  real<lower=0> sigma_a;
  real mu_b;
  real<lower=0> sigma_b;
}
model {
  // priors
  mu_a ~ normal(0,100);
  mu_b ~ normal(0,10);
  sigma_a ~ cauchy(0,1);
  sigma_b ~ cauchy(0,1);
  for (j in 1:n_group){
    a[j] ~ normal(mu_a,sigma_a);
    b[j] ~ normal(mu_b,sigma_b);
  }
  sigma ~ normal(0,100);
  // likelihood
  for(i in 1:n){
    y[i] ~ normal( a[ group[i] ] + b[ group[i] ] *x[i] , sigma); 
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

fit_partial  = sampling(stan_model_partial,
                        data=data,
                        chains=3,
                        iter=5000,
                        warmup=2000,
                        control=list(adapt_delta=0.99)
)

print(fit_partial, digits=3, probs=c(0.025, 0.975))

plot(fit_partial)
# plot(fit_partial, pars="a")
plot(fit_partial, pars="a")
plot(fit_partial, pars="b")
plot(As.mcmc.list(fit_partial)) # from coda package

# pairs(fit_partial, pars=c("b","mu_b","sigma_b"))

posterior=as.matrix(fit_partial)

#------------------------------------------------------------------------------
# predictions / credible intervals
#------------------------------------------------------------------------------

# Again, we can generate predictions and compute credible intervals (for the deterministic part of the model)
# (Here: 90% credible intervals)

x.pred = seq(from=min(df$body_mass_g_norm), to=max(df$body_mass_g_norm), by=0.1)

par(mfrow=c(3,1))

for (i in 1:3){
  df.sub=subset(df, df$group==i)
  plot(df.sub$body_mass_g_norm,  df.sub$bill_length_mm,
       xlim=range(df$body_mass_g_norm),
       ylim=range(df$bill_length_mm),
       main=levels(df$species)[i]
  )
  
  y.cred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))
  for(j in 1:nrow(posterior)){
    y.cred[j, ] = posterior[j,paste0("a[",i,"]")] + posterior[j,paste0("b[",i,"]")]*x.pred 
  }
  
  y.cred.mean = apply(y.cred, 2, function(x) mean(x)) 
  lines(x.pred, y.cred.mean, col="red", lwd=2)
  
  y.cred.q05 = apply(y.cred, 2, function(x) quantile(x, probs=0.05)) 
  lines(x.pred, y.cred.q05, col="red", lwd=2, lty=2)
  
  y.cred.q95 = apply(y.cred, 2, function(x) quantile(x, probs=0.95)) 
  lines(x.pred, y.cred.q95, col="red", lwd=2, lty=2)
}

#------------------------------------------------------------------------------
# Further reading
#------------------------------------------------------------------------------

# read about complete pooling, partial pooling, no pooling and shrinkage:

# https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/
  
  
