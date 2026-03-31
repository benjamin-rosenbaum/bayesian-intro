# Practical 5: Generalized linear models
# by Benjamin Rosenbaum

# We learn how to use some classic GLMs, and a distributional model. We focus on 
# prior specifications and appropriate posterior predictive checks, including 
# the DHARMa package. Sometimes we have to apply ggplot tricks for good conditional 
# effects plots. 

# General disclaimer: Never use LOOIC / WAIC / AIC / BIC / DIC for comparing models 
# for continuous responses (Normal, Lognormal, Gamma, Beta, ...) to models for 
# discrete responses (Binomial, Poisson, Negative Binomial, ...)
# Continuous distributions are based on probability density functions, while 
# discrete distributions have probability mass functions, which you can't compare.

rm(list=ls())
library("brms")
library("ggplot2")
library("arm")
library("emmeans")
library("performance")
library("Data4Ecologists")
library("DHARMa.helpers")

setwd("~/Nextcloud/Teaching brms 2025 2/Practical_05")

# Count variables ----------------------------------------------------------------------

## Poisson

# We start with a simple comparison of 2 group means (frequentists would use a 
# t-test), but the response are counts (number of slugs found per plot, 2 field 
# types, from Data4Ecologists package).

data(slugs)
ggplot(slugs, aes(field, slugs)) + geom_jitter(alpha=0.5, width=0.1, height=0.05)

# deterministic part: mu[field[i]] --> just mu1 or mu2
# stochastic part:    y[i] ~ Poisson(mu[field[i]])
# don't use a log-link for this simple model, just "identity link"

default_prior(slugs ~ field,
              family = poisson(link=identity),
              data = slugs)

# just provide a vague prior for effect (difference in means, dummy-coded)
# and ensure that intercept is positive with lower boundary (`lb=0`)

fit.slugs = brm(slugs ~ field,
                family = poisson(link=identity),
                prior = 
                  prior(normal(0,2), class=b) +
                  prior(student_t(3, 1, 2.5), class=Intercept, lb=0),
                data = slugs)

# check convergence

plot(fit.slugs)
summary(fit.slugs, prior=TRUE)

# posterior predictions

plot(conditional_effects(fit.slugs))
plot(conditional_effects(fit.slugs), 
     points=T, 
     point_args=list(alpha=0.3, width=0.1, height=0.05)
)

dh_check_brms(fit.slugs) # DHARMa

# DHARMa checks don't look so good. Poisson too strict because it assumes 
# mean = variance. 

## Negative Binomial = overdispersed Poisson

# An additional precision parameter called `shape` controls dispersion: 
# For large values, we get the Poisson, for small values the Negative Binomial 
# is an overdispersed Poisson.

# It is classically specified by `family=negbinomial()`, but since the identity 
# link is not a standard choice we have to "convince" brms by using 
# `family=brmsfamily("negbinomial, )"`.

default_prior(slugs ~ field,
              family = brmsfamily("negbinomial", link="identity"), 
              data = slugs)

fit.slugs.nb = brm(slugs ~ field,
                   family = brmsfamily("negbinomial", link="identity"),
                   prior = 
                     prior(normal(0, 2), class=b) +
                     prior(student_t(3, 1, 2.5), class=Intercept, lb=0),
                   data = slugs)

# Check convergence

summary(fit.slugs.nb, prior=TRUE)
plot(fit.slugs.nb)

# shape is small, so the data is really overdispersed. 

# Posterior predictions

plot(conditional_effects(fit.slugs.nb), 
     points=T, 
     point_args=list(alpha=0.3, width=0.1, height=0.05)
)
dh_check_brms(fit.slugs.nb) # DHARMa

# The distributional assumptions are OK now. 

LOO(fit.slugs, fit.slugs.nb)

# Model comparison also tells us that the Negative Binomial is more appropriate.

fixef(fit.slugs) |> round(3)
fixef(fit.slugs.nb) |> round(3)

# The effect size (difference in group-means) has a larger standard deviation 
# in the Negative Binomial. In different words: we would have been overly confident 
# about the effect size in the Poisson model.

hypothesis(fit.slugs, "fieldRookery>0")
hypothesis(fit.slugs.nb, "fieldRookery>0")


# Binomial ----------------------------------------------------------------------

# data from Quian (2016) Environmental and Ecological Statistics with R

# mice infected with parasites of the genus Cryptosporidium 
# fit a dose-response model
# N = number of inoculated mice
# Y = number of infected mice
# Dose = number of parasite oocysts used for inoculation
# Which minimal dose is needed to infect at least 75% of mice?

# deterministic part: logit(p) = a + b*Dose
# stochastic part:    Y[i] ~ Binomial(N[i],p[i])

df = read.csv("https://raw.githubusercontent.com/songsqian/eesR/refs/heads/master/R/Data/cryptoDATA.csv",
              sep=" ")
df.sub = subset(df, Source=="Finch")
df.sub$Y = as.integer(df.sub$Y)

ggplot(df.sub, aes(Dose, Y/N)) + geom_point(alpha=0.5)

# can also use scale() as predictor. 
# priors and parameters on z-scale, plots / conditional_effects on original scale

default_prior(Y | trials(N) ~ scale(Dose),
              family = binomial(link=logit),
              data = df.sub )

# add a weakly informative prior on (linear scale) slope. effect is expected to be pos

fit.crypto = brm(Y | trials(N) ~ scale(Dose),
                 family = binomial(link=logit),
                 prior = prior(normal(1,1), class=b),
                 data = df.sub )

# check convergence

summary(fit.crypto, prior=T)
plot(fit.crypto)

# conditional_effects plots predictions for N=1 only.
# "points=T" just adds Y values, but we need Y/N values

plot(conditional_effects(fit.crypto, effects="Dose"), points=T)

# here we go:

p1 = plot(conditional_effects(fit.crypto, effects="Dose"), 
          points=FALSE, plot=FALSE)
p1[[1]] + geom_point(data=df.sub, 
                     aes(Dose, Y/N), 
                     alpha=0.6, size=1.5,
                     inherit.aes=FALSE)

# posterior predictions

pp_check(fit.crypto, ndraws=100) 
pp_check(fit.crypto, type="scatter_avg")

# not very pretty, but we didn't measure any other predictors to improve
# could try beta_binomial instead of binomial

fitted = fitted(fit.crypto)
residuals = residuals(fit.crypto)
binnedplot(fitted, residuals)
dh_check_brms(fit.crypto)

# What parasite dose is needed to get at least 75% of mice infected?

p1[[1]] + geom_hline(yintercept=0.75, linetype=2, col="red")

# for Dose=191, on average 75% are infected
# but we are just 50% sure that more than 75% are infected 

fitted(fit.crypto, 
       newdata = data.frame(Dose=191, N=1))
fitted.191 = fitted(fit.crypto, 
                    newdata = data.frame(Dose=191, N=1),
                    summary=FALSE)
hist(fitted.191)
abline(v=0.75, col="red", lwd=2)

# for a Dose=229, the lower 2.5% quantile of fitted values exceeds 75%
# --> we are 97.5% sure that more than 75% are infected

fitted(fit.crypto, 
       newdata=data.frame(Dose=229, N=1), 
       probs=0.025)
fitted.229 = fitted(fit.crypto, 
                    newdata = data.frame(Dose=229, N=1),
                    summary=FALSE)
hist(fitted.229)
abline(v=0.75, col="red", lwd=2)


# Distributional model ---------------------------------------------------------

# daily growth of young chicken (from datasets package). 
# only early exponential growth phase here, first 2 weeks
# (otherwise nonlinear "von Bertalanffy growth function")

# Question: What is the daily growth rate?
  
data("ChickWeight")

ggplot(ChickWeight, aes(Time, weight)) + 
  geom_jitter(alpha=0.5, width=0.1, height=0.05)

ggplot(ChickWeight, aes(Time, log(weight))) + 
  geom_jitter(alpha=0.5, width=0.1, height=0.05)

ChickWeight = subset(ChickWeight, Time<15)

# simple linear regression for log(weight). no link function, just transformed response

# deterministic part: mu = a + b*time
# stochastic part:    log.weight ~ normal(mu, sigma) 

fit.growth.1 = brm(log(weight) ~ Time,
                   prior = prior(normal(0,1), class=b),
                   data = ChickWeight)

summary(fit.growth.1)
plot(fit.growth.1)

# posterior predictions (for deterministic part mu)

plot(conditional_effects(fit.growth.1), 
     points=T, 
     point_args=list(alpha=0.4))

# posterior predictions (predicted data, includes stochastic part)

plot(conditional_effects(fit.growth.1, method="posterior_predict"), 
     points=T, 
     point_args=list(alpha=0.4))
# --> prediction intervals don't represent variance of data well

pp_check(fit.growth.1, type="scatter_avg")
pp_check(fit.growth.1, ndraws=100)
check_model(fit.growth.1, check=c("linearity","homogeneity","qq","normality"))
dh_check_brms(fit.growth.1) # DHARMa
# --> linear model assumptions not satisfied

# next: distributional model, 
# not only mu, but also sigma is a function of the predictor time.

# deterministic part: mu = a + b*time
#                     log(sigma) = c + d*time (log-link to keep sigma positive)
# stochastic part:    log.weight ~ normal(mu, sigma) 

default_prior(bf(log(weight) ~ Time,
                 sigma ~ Time),
              data = ChickWeight)

fit.growth.2 = brm(bf(log(weight) ~ Time,
                      sigma ~ Time),
                   prior = prior(normal(0,1), class=b, coef=Time) +
                           prior(normal(0,1), class=b, dpar=sigma),
                   data = ChickWeight)

summary(fit.growth.2, prior=T)
plot(fit.growth.2)

# posterior predictions (for deterministic part mu)

plot(conditional_effects(fit.growth.2), 
     points=T, 
     point_args=list(alpha=0.4))

# posterior predictions (predicted data, includes stochastic part)

plot(conditional_effects(fit.growth.2, method="posterior_predict"), 
     points=T, 
     point_args=list(alpha=0.4))
# --> looks much better

pp_check(fit.growth.2, type="scatter_avg")
pp_check(fit.growth.2, ndraws=100)
# We don't use check_model() here, because it's not a linear model anymore.

dh_check_brms(fit.growth.2) 
# --> DHARMa checks OK-ish

LOO(fit.growth.1, fit.growth.2)
# --> distributional model represents data much better

# Question: What's the average daily growth rate?
# log.weight=a+b*time --> weight=exp(a)*exp(b)^time
# exp(b) is daily growth rate. use posterior distribution of exp(b) to compute mean

post = as_draws_df(fit.growth.2)
post$rate = exp(post$b_Time)
hist(post$rate)

mean(post$rate)
quantile(post$rate, prob=c(0.05, 0.95))

# log link on sigma kept it positive, but cause exponential increase with time
# alternative model without sigma log-link --> sigma increases linearly with time
# identity link, but sigma's intercept and slope positive

# fit.growth.3 = brm(bf(log(weight) ~ Time,
#                       sigma ~ Time),
#                    family = brmsfamily("gaussian", link_sigma = "identity"),
#                    prior = 
#                      prior(normal(0,1), class=b, coef=Time) +
#                      prior(normal(0,1), class=Intercept, dpar=sigma, lb=0) +
#                      prior(normal(0,1), class=b, dpar=sigma, lb=0),
#                    data = ChickWeight)

# nonlinear model --------------------------------------------------------------

# TO DO !!!


# exercise: logistic -----------------------------------------------------------

# From Qian, S. (2016) Environmental and Ecological Statistics with R 

# Seed predation by rodents
# Seeds were placed in gauze bags in 4 topographic locations (topo)
# Visited in 6 sampling campaigns (time)
# Q: Does predation rate reach a stable value in the last three time periods?

# Start with Predation ~ log(seed.weight)
# Test if there is a difference between locations (+ topo)
# Also use time as factorial predictor (+ time)

# Question: Does predation rate reach a stable value in the last three time periods?
  
df = read.csv("https://raw.githubusercontent.com/songsqian/eesR/refs/heads/master/R/Data/seedbank.csv",
              sep=",")

# data preparation
df = subset(df, Predation %in% 0:1)
df$Predation = as.integer(df$Predation)
df$time = as.factor(df$time)
df$topo = as.factor(df$topo)
df$seed = scale(log(df$seed.weight))
head(df)

ggplot(df, aes(log(seed.weight), Predation)) + 
  geom_jitter(height=0.05, width=0.05, alpha=0.2)

ggplot(df, aes(time, Predation)) + 
  geom_jitter(height=0.05, width=0.05, alpha=0.2)

