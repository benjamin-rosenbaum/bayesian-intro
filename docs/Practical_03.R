# Practical 3: Priors & posteriors
# by Benjamin Rosenbaum

# We learn about the posterior distribution and posterior predictions.
# Posterior predictive checks are used for model evaluation.
# Prior predictive checks are used to see if prior makes sense.

rm(list=ls())
library("brms")
library("bayesplot") 
library("performance")
library("ggplot2")
library("ecostats")
try(dev.off())

setwd("~/Nextcloud/Teaching brms/Practical_03")

# ecostats::globalPlants -------------------------------------------------------

# This dataset contains observations of plant height vs latitude.
# There is a negative relationship expected.
# We perform a regression of log height vs latitude.

# Deterministic part: mu = a+b*lat
# Stochastic part:    logheight ~ Normal(mu,sigma)

data(globalPlants)
str(globalPlants)
plot(globalPlants$lat, globalPlants$height)
plot(globalPlants$lat, log(globalPlants$height))

fit1 = brm(log(height) ~ lat,
           prior = prior(normal(-0.05,0.02), class=b, coef=lat),
           data=globalPlants )

# check  convergence

summary(fit1, prior=TRUE)
plot(fit1)
plot(conditional_effects(fit1), points=T)

# ++ posterior distribution -------------------------------------------------------

# What exactly is the posterior? It's just a collection of (multivariate) samples 
# of model parameters. They can be extracted as a matrix or dataframe
# (rows=samples, columns=parameters)

post = as_draws_matrix(fit1) 
str(post)
head(post)

post = as_draws_df(fit1)
str(post)
head(post)

# For example, we can look at the distribution of the intercept, which is the 
# fitted log(height) at the equator (lat=0).

mcmc_hist(fit1, pars="b_Intercept") 
hist(post$b_Intercept) 

# The posterior contains 2 variables for intercept? What's going on?
# The other variable `Intercept` is used by brms internally for model fitting 
# and describes the intercept for mean centered predictors. It is the fitted 
# log(height), where the predictor is at its empirical mean (lat=mean(lat))

mcmc_hist(fit1, pars="Intercept") 
hist(post$Intercept) 

# posterior is multivariate!
# classical neg. correlation of int and slope. but this is not a huge problem

pairs(fit1, variable=c("b_Intercept", "b_lat", "sigma"))
pairs(fit1, variable=c("b_Intercept", "b_lat", "sigma"), off_diag_fun="hex")

# correlation disappears if parametrized with mean-centered predictor
# makes MCMC more efficient, but model is equivalent

pairs(fit1, variable=c("Intercept", "b_lat", "sigma"), off_diag_fun="hex")

# ++ posterior predictions ------------------------------------------------------

# For any given value of the predictor, predictions can be computed from the 
# samples of the posterior. If we just use the deterministic model part, this is 
# called the fitted distribution, when we additionally include the stochastic 
# model part, this is called the predictive distribution.

# fitted distribution (deterministic part)
# posterior predictions for mu=a+bx 
# each sample (row in the posterior) is a regression line

head(post[, 1:3])
plot(conditional_effects(fit1), points=TRUE)
plot(conditional_effects(fit1, prob=0.80), points=T)
plot(conditional_effects(fit1, prob=0.99), points=T)
plot(conditional_effects(fit1, spaghetti=TRUE, ndraws=100), points=TRUE)
# 95% credible intervals display uncertainty of mu

# predicted distribution (deterministic + stochastic part)
# posterior predictions for data y=mu+eps, eps~Normal(0,sigma)

plot(conditional_effects(fit1, method="posterior_predict"), points=TRUE)
# 95% prediction intervals should contain ~95% of datapoints

# can extract credible & prediction intervals for all observations (datapoints)

fitted(fit1) |> round(3)
predict(fit1) |> round(3)

# or make predictions for new predictor value

fitted(fit1, newdata=data.frame(lat=60)) |> round(3)
predict(fit1, newdata=data.frame(lat=60)) |> round(3)

# inference with predictions
# use hypothesis function for inference on parameters or predictions

hypothesis(fit1, "lat<0")

# but what about: whats the expected difference in logheight between 40 and 60 lat?

fitted(fit1, newdata=data.frame(lat=40)) |> round(3)
fitted(fit1, newdata=data.frame(lat=60)) |> round(3)

# we need posterior distribution of mu(60)-mu(40)

mu40 = posterior_epred(fit1, newdata=data.frame(lat=40))
mu60 = posterior_epred(fit1, newdata=data.frame(lat=60))

cbind(mu40,mu60) |> head()

hist(mu40)
hist(mu60)
hist(mu40-mu60)
mean(mu40-mu60)
sd(mu40-mu60)
quantile(mu40-mu60, probs=c(0.05, 0.90))
sum(mu40>mu60)/length(mu40)

# alternatively use hypothesis function (a bit complicated, but possible)

mus = fitted(fit1, 
             newdata=data.frame(lat=c(40,60)),
             summary=F)
mus = as.data.frame(mus)
names(mus)=c("mu40","mu60")
hypothesis(mus, "mu40>mu60")
plot(hypothesis(mus, "mu40>mu60"))

# ++ posterior predictive checks -----------------------------------------------

# How well does the model fit the data?

# bayes R2 (is calculated a bit different from freq R2, but means the same)

bayes_R2(fit1)

# use predictive distribution to see if distribution of response values is replicated

pp_check(fit1, ndraws=100)

# classical observed vs predicted plot

pp_check(fit1, type="scatter_avg")

# some analyses from the performance package for linear models
# (more on that in next session)

check_model(fit1, check=c("linearity","homogeneity","qq","normality"))

# ++ prior predictive checks ---------------------------------------------------

# Now that we understand how to do posterior predictive checks, we can also do 
# prior predictive checks.

# Unless we have prior information from previous studies / experiments, how do 
# we know if we have a good / meaningful prior? Parameters sampled from the prior
# distribution should generate predictions which are roughly in the range of 
# observed data. Sure, the information contained in the data goes into the model 
# via the likelihood, but priors that generate predictions of a total different 
# magnitude as the data do not make much sense. 

get_prior(fit1)

# instead of sampling from the posterior ~ prior * likelihood, 
# we just sample from the prior. 
# why? check if predictions from prior are approximately in line with data
# for lin models this may be obvious, but for GLMs etc more important & helpful

fit1.prior = brm(log(height) ~ lat,
                 prior = prior(normal(-0.05,0.02), class=b, coef=lat),
                 sample_prior = "only",
                 data=globalPlants )

summary(fit1.prior, prior=TRUE)

plot(conditional_effects(fit1.prior), points=TRUE)
plot(conditional_effects(fit1.prior, spaghetti=TRUE, ndraws=100), points=TRUE)
pp_check(fit1.prior, ndraws=100)

# --> prior predictions covering the data well. same order of magnitude
# --> this is a vaguely informative prior

# could also tighten the prior a bit (especially brms default for intercept),
# but not that important here

fit2.prior = brm(log(height) ~ lat,
                 prior = c(prior(normal(-0.05,0.02), class=b, coef=lat),
                           prior(student_t(3, 1.1, 1), class=Intercept)),
                 sample_prior = "only",
                 data=globalPlants )

summary(fit2.prior, prior=TRUE)

plot(conditional_effects(fit2.prior), points=TRUE)
plot(conditional_effects(fit2.prior, spaghetti=TRUE, ndraws=100), points=TRUE)
pp_check(fit2.prior, ndraws=100)

# --> a bit more in line with the data

# fit model with the new prior

fit2 = brm(log(height) ~ lat,
           prior = c(prior(normal(-0.05,0.02), class=b, coef=lat),
                     prior(student_t(3, 1.1, 1), class=Intercept)),
           data=globalPlants )

# here, prior choice affects the posterior only marginally

fixef(fit1)
fixef(fit2)

# ecostats::waterQuality -------------------------------------------------------

# Test if river health deteriorates downstream (increased catchment area)

# Analyze the data of water quality in several lakes vs area (log catchment area)
# Now perform previous steps in the correct logical order 
# (1) check /plot data. scale predictor --> scale()
# (2) assign priors
# (3) fit model
# (4) check convergence
# (5) check model fit
# (6) inference: 
#     mean diff in water quality between mean and mean+1sd of pred, and its 90%-CI

# (1) check data ---------------------------------------------------------------

data("waterQuality")
str(waterQuality)
plot(scale(waterQuality$logCatchment), waterQuality$quality)
waterQuality$log.catch.z = scale(waterQuality$logCatchment)
