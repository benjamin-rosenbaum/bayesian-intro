rm(list=ls())
library(rstan)
library(coda)
library(BayesianTools)
# setwd("~/Desktop/teaching Bayes")

# Aim: learn what a random effect is.
# Learn about hierarchical levels of parameters. Bayesian statistics is great for hierearchical models!
# Think in terms of "no pooling", "partial pooling" and "complete pooling" of parameters
# rather than "random effects" and "fixed effects".

rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

#------------------------------------------------------------------------------
# load data 
#------------------------------------------------------------------------------

# We load an example dataset from Kruschke "Doing Bayesian Data Analysis".
# It contains information on the lifespan (Longevity) of male fruitflies depending on their sexual activity.
# Sexual activity was experimentally manipulated by different numbers of pregnant or virgin female fruitflies (CompanionNumber).

df = read.csv("~/git/bayesian-intro/data/FruitflyDataReduced.csv")
head(df)
str(df)

boxplot(Longevity ~ CompanionNumber, data=df,
        ylab="Longevity",
        xlab="Companions ID",
        col="grey")

plot(0, 0, xlim =c(1,5), ylim = range(df$Longevity), type = "n",
     ylab="Longevity",
     xlab="Companions ID")
points(Longevity ~ jitter(as.numeric(CompanionNumber), factor=0.2), data=df)

# The following approach is similar to frequentist ANOVA.

#------------------------------------------------------------------------------
# no pooling / fixed effects model 
#------------------------------------------------------------------------------

# We want to estimate the mean longevity per group, and test if there are differences between groups.
# The statistical model is identical to yesterday's "t-test" model, but now there are more than 2 groups:
# y_i ~ normal(mu[group_i], sigma), i=1,...,n

# Note that there are no assumptions on mu, i.e. the means of the different groups are estimated independently.
# "group" is treated as a fixed effect. 
# Because there can be misunderstandings when using the terms "fixed effect" / "random effects", 
# this is also called "no pooling", i.e. no information on the mu_j is pooled across the groups. Keep that in mind for later!

# Here, we used the same variance for each group (sigma), but we could also use different sigmas.

# When preparing the data for Stan, note that we use "as.integer()" to code the factorial variable "group",
# so we can use it as an index.

data = list(y = df$Longevity,
            group = as.integer(df$CompanionNumber),
            n = nrow(df),
            n_group = 5)

data

# We use (very) weakly informative prior information on the mu_j, but note that the priors for mu_j are assigned independently!

stan_code_nopool = '
data {
  int n;
  int n_group;
  real y[n];
  int group[n];
}
parameters {
  real<lower=0> mu[n_group];  
  real<lower=0> sigma;  
}
model {
  // priors
  for (j in 1:n_group){
    mu[j] ~ normal(0,100);
  }
  sigma ~ normal(0,100);
  // likelihood
  for(i in 1:n){
    y[i] ~ normal(mu[group[i]], sigma);
  }
}
'

stan_model_nopool = stan_model(model_code=stan_code_nopool)
# save(stan_model_nopool, file="stan_code_nopool.RData")
# load("stan_code_nopool.RData")

fit_nopool  = sampling(stan_model_nopool,
                       data=data,
                       chains=3,
                       iter=2000,
                       warmup=1000
)

print(fit_nopool, digits=3, probs=c(0.025, 0.975))

plot(fit_nopool)
# plot(fit_nopool, pars="mu")
# plot(As.mcmc.list(fit_nopool)) # from coda package

# As in the "t-test" example, we can look at the individual differences between groups, also called "contrasts".
# E.g., the posterior distribution of mu4-mu5

posterior_nopool = as.matrix(fit_nopool)

contrast = posterior_nopool[,"mu[4]"]-posterior_nopool[,"mu[5]"]

plot(density(contrast))
abline(v=0, col="red", lwd=2)

#------------------------------------------------------------------------------
# partial pooling / random effects model 
#------------------------------------------------------------------------------

# Until now, we assumed that the different group means are completely independent.
# However, we can assume that the longevity of fruitflies is not completely independent.
# We can assume there is a joint mean and each group mean differs from this joint mean with a group-level standard deviation.
# And we can explicitely model that:

# y_i ~ normal(mu[group_i], sigma), i=1,...,n (n observations)
# mu_j ~ normal(mu_total, sigma_total), j=1,...,m (m groups)

# sigma describes within-groups variation of response values,
# sigma_total describes between-groups variation of group means.

# Here, group is treated as a random effect. 
# Because some information is shared / pooled across groups, this is also called "partial pooling".

# mu_group and sigma_group are also parameters that are estimated from the data.
# Because mu_j's distribution depends on mu_total and sigma_total, this is a hierarchical model. 

# Like all other parameters, they will be assigned a prior distribution, too.
# For the joint mean mu_total, we can use the same expectation we used for the group means before.
# For the between-groups standard deviation sigma_total, it is standard procedure to use half-Chauchy distributions
# (more heavy-tailed than normal distribution)

# The partial pooling model differs from the no pooling model above only in the joint distribution of the mu_j 
# and additional parameters mu_total and sigma_total.

stan_code_partpool = '
data {
  int n;
  int n_group;
  real y[n];
  int group[n];
}
parameters {
  real<lower=0> mu[n_group];  
  real<lower=0> sigma;  
  real<lower=0> mu_total;
  real<lower=0> sigma_total;
}
model {
  // priors
  mu_total ~ normal(0, 100);
  sigma_total ~ cauchy(0, 10);

  for (j in 1:n_group){
    mu[j] ~ normal(mu_total, sigma_total);
  }
  sigma ~ normal(0,100);
  
  // likelihood
  for(i in 1:n){
    y[i] ~ normal(mu[group[i]], sigma);
  }
}
'

stan_model_partpool = stan_model(model_code=stan_code_partpool)
# save(stan_model_partpool, file="stan_code_partpool.RData")
# load("stan_code_partpool.RData")

fit_partpool  = sampling(stan_model_partpool,
                         data=data,
                         chains=3,
                         iter=2000,
                         warmup=1000
)

print(fit_partpool, digits=3, probs=c(0.025, 0.975))

plot(fit_partpool)
# plot(fit_partpool, pars="mu")
# plot(As.mcmc.list(fit_partpool)) # from coda package

#------------------------------------------------------------------------------
# Comparison
#------------------------------------------------------------------------------

# Now we compare the results of the "no pooling" and the "partial pooling" models.
# The means and 95% confidence intervals of the mu_j can directly be extracted in a summary.

summary_nopool = summary(fit_nopool)$summary
summary_nopool

summary_partpool = summary(fit_partpool)$summary
summary_partpool

plot(0, 0, xlim =c(1,5), ylim = range(df$Longevity), type = "n",
     ylab="Longevity",
     xlab="Companions ID")
points(Longevity ~ jitter(as.numeric(CompanionNumber), factor=0.2), data=df,
       col="grey")

points(summary_nopool[1:5, "mean"], pch="+", col="blue", cex=2)
points(summary_nopool[1:5, "2.5%"], pch="-", col="blue", cex=2)
points(summary_nopool[1:5, "97.5%"], pch="-", col="blue", cex=2)

points(summary_partpool[1:5, "mean"], pch="+", col="red", cex=2)
points(summary_partpool[1:5, "2.5%"], pch="-", col="red", cex=2)
points(summary_partpool[1:5, "97.5%"], pch="-", col="red", cex=2)

abline(h=mean(df$Longevity), col="grey")

legend("topright", legend=c("no pooling","partial pooling"), pch=c("+","+"), col=c("red","blue"), bty="n")

# The difference between "no pooling" and "partial pooling" estimates is that 
# extreme values tend to be pulled towards the joint mean by partial pooling (shrinkage).
# Statistical power is borrowed across groups.
# This can be helpful if some groups have a low number of observations.

# In addition to both models, there is also "complete pooling". Here, all information across groups would be pooled:
# y_i ~ normal(mu, sigma)
# There is no effect of group included, i.e. all groups share the same mean. 
# "partial pooling" is a compromise between "no pooling" and "complete pooling".


