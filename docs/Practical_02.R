# Practical 2: The brms package
# by Benjamin Rosenbaum

# We learn about the brms package and how to fit simple regression models.
# Focus on: model output, convergence checks. 
# Also some basic infos on priors and model predictions. 

rm(list=ls())
library("brms")
library("ggplot2")
library("bayesplot") 
library("loo") 
try(dev.off())

setwd("~/Nextcloud/Teaching brms/Practical_02")

# (1) linear regression --------------------------------------------------------

# We start with our deer population and the simple weight~age example

# Q: What's the average growth per year? (Slope in age)
# Deterministic part: mu = a+b*age
# Stochastic part:    weight ~ Normal(mu,sigma)

data = data.frame(weight = c(104, 120, 118, 115, 99, 110, 102),
                  age    = c(10, 12, 11, 11, 9, 11, 10))
plot(data$age, data$weight)

# basic brms functions----------------------------------------------------------

# Instead of "lm()", we use the "brm()" function. The formula notation is designed 
# to be identical to lm, glm, lme4 (with few exceptions) 

fit1 = brm(weight ~ age, data=data)

# Looking at the summary table, we get a lot of infos: 
# brm by default uses 4 chains, each with 1000 warmup and 1000 sampling iterations.
# The first thing you should look at are not parameter estimates, but Rhat and ESS.
# These indicate if the MCMC converged and the posterior distribution is properly
# sampled. Check if Rhat<1.1 and compare ESS to total number of draws.

summary(fit1)

# Additionally, you should do a visual inspection of the MCMC. You get a histogram
# and a traceplot per parameter, which should look like a fuzzy caterpillar

plot(fit1)

# You can change the color palette if you like

color_scheme_set("viridisA")
plot(fit1)

# You can also specify to display just some selected parameters.
# Parameters of the deterministic model part begin with "b_", 
# the residual standard deviation is "sigma"

color_scheme_set("blue")
plot(fit1, variable=c("b_Intercept", "b_age"))
plot(fit1, variable=c("sigma"))

# Histograms and traceplots can be plotted individually

mcmc_trace(fit1, pars=c("b_Intercept", "b_age"))
mcmc_hist(fit1, pars=c("b_Intercept", "b_age"))

# All the brms plots are done withh ggplot2, so you can extract & modify them

plot1 = mcmc_hist(fit1, pars=c("b_age"))
plot1 + xlim(-10,30) 

# Now that we have verified MCMC convergence, let's get back to the model itself.
# If not specified otherwise, brms uses a normal distribution for the residuals:
# family=gaussian(). We get some information on estimated parameters (mean etc).

summary(fit1)

# In this simple 1-predictor regression, model predictions are easily plotted vs data.
# `conditional_effects()` is a powerful function which we will use throughout the course.

plot(conditional_effects(fit1))
plot(conditional_effects(fit1), 
     points=TRUE)

# Again, this generates a ggplot object which can be modified, with some options
# in the plot function, or full ggplot options if you save the object

plot(conditional_effects(fit1), 
     points=TRUE, point_args=list(col="red", alpha=0.5))

# Note that we here only plot the uncertainty of the deterministic model part mu,
# more on that tomorrow. "fitted" computes predictions of mu for each datapoint.

fitted(fit1)

# The brms package does not only offer model fitting via MCMC, it also has a lot
# of functions for model analysis and is compatible with a lot of other packages 
# (e.g. emmeans). We will learn about some of these in the next days.

methods(class="brmsfit")


# brms specifications & priors -------------------------------------------------

#  When we compare the results to frequentist lm-model, the slope is pretty close
#  but there's ~0.5 difference in intercepts. So why are they different? 
#  What about priors, did we specify any?

fit1.lm = lm(weight ~ age, data=data)
summary(fit1.lm)
summary(fit1)

# The brm function has A TON OF specifications, which we did not specify in the 
# simple brm(weight~age, data=data) model. So brms uses default values.

?brm

# E.g. we can specify the number of chains & iterations manually.
# Per default, half of the iterations are used for warmup and are discared from 
# the posterior sample.

fit2 = brm(weight ~ age,
           data = data,
           chains = 4,
           iter = 5000
           )

# With a larger number of samples, we expect a more accurate approximation of the
# true posterior. Parameter means usually are quite correct even for low numbers,
# while outer quantiles (e.g. 90%, 95%) require larger numbers of samples

summary(fit2)
plot(fit2)

# OK, but what about priors? Are there any defaults used? 
# We can check the defaults for any model with `default_prior()`. 
# The model does not have to be fitted, just model formula and data must be specified.

default_prior(weight ~ age,
              data = data)

# Alternatively, you can display the priors of any fitted model.

prior_summary(fit2)

# Since we did not specify any priors, both displays are the same here.

# This table can be a bit confusing, but look at the column "class":
# "b" is for effects / slopes. The first line tells you if there is a prior used
# for ALL effects, which is not the case (prior=flat). Second line is the prior for
# a specific coefficient (coef=age), there's also no prior specified.
# But brms chooses a prior for "Intercept" and for the residual sdev "sigma". 
# These are automatically generated from the mean and the spread of the response.
# Note that internally, the brms machine uses mean-centered predictors.
# The "Intercept" parameter (and its prior) are based on mean-centered variables.
# What's displayed in the model summary is actually "b_Intercept" which is the 
# intercept parameter transformed to the original, non mean-centered scale.

# A short form is presented in the summary for prior=TRUE

summary(fit2, prior=TRUE)


# Unless necessary, I would leave the brms defaults for Intercept & sigma.
# However, you should choose a prior for the slope, which currently has none.

# This would set a prior for all slopes (if you have >1 predictors)

my_priors = prior(normal(5,1), class=b)

# For setting a prior for a specific predictor, you specify it in "coef". 
# Since this model only has 1 predictor, both formulations are the same.

my_priors = prior(normal(5,1), class=b, coef=age)

fit3 = brm(weight ~ age,
           prior = my_priors,
           data = data,
           chains = 4,
           iter = 5000
)

# same as:
# fit3 = brm(weight ~ age,
#            prior = prior(normal(5,1), class=b, coef=age),
#            data = data,
#            chains = 4,
#            iter = 5000
# )

summary(fit3, prior=TRUE)
plot(fit3)

plot(conditional_effects(fit3), 
     points=TRUE)

# lin.reg., model analysis -----------------------------------------------------

# Only after we checked MCMC convergence, we can go to the next step:
# Model evaluation / model checking. How well does our model describe the data?

# A classical visual tool is observed vs predicted, which also works if you have 
# multiple predictors. 

pp_check(fit3, "scatter_avg")

# bayes_R2 is the amount of explained variation. It's computation is a bit 
# different from the classical freq. R2, but conceptually it means the same.

bayes_R2(fit3)

# More on that tomorrow, e.g. checking model assumptions.

# lin.reg., inference ----------------------------------------------------------

# OK, so we know that (a) MCMC converged and (b) model describes the data well.
# Only now can we make inference, i.e. quantitative statements about research questions.
# The summary already tells us mean and 95-CIs for the slope (growth per year of age)

summary(fit3)

# Different Credible intervals can be chosen in the summary, e.g. 80%-CI

summary(fit3, prob=0.80)

# Or you can extract specific quantiles of parameter estimates

fixef(fit3)
fixef(fit3, probs=c(0.25, 0.5, 0.75))

# In a frequentist analysis you would want to know if the effect of age is "significant":
# p-values quantify the probability of observing such a pattern the data 
# if the null hypothesis (b_age=0) was true. (p small -> reject H0)
# Here, we can just calculate the probability that the slope is positive
# "Post.Prob" is the value of interest. It's =1 because all samples of slope were positive

hypothesis(fit3, "age>0")
plot(hypothesis(fit3, "age>0"))

# You can test all kinds of hypotheses for the parameters! If we were interested
# in the question if growth per year is bigger than 4, just put it in the hypothesis:

hypothesis(fit3, "age>4")
plot(hypothesis(fit3, "age>4"))

# It's transformed in the equivalent formulation "age-4>0", and this is the
# probability distribution which is actually plotted. 
# The posterior probability is 0.96, which is also the integral 
# (area under the curve) for "age-4>0"

plot1 = plot(hypothesis(fit3, "age>4"), plot=FALSE)
plot1[[1]] + geom_vline(xintercept=0) 


# (2) survival rate ----------------------------------------------------------------

data = data.frame(total = c(22,22,29,21,25,30,24,23,25,28),
                  survived = c(19,14,23,19,20,18,15,16,18,15))

# Exercise: Analyze the model for 
# deterministic: mu=theta (in [0,1])
# stochastic:    yi ~ Binomial(N_i,theta)
# (just fitting an intercept)

# brm(survived | trials(total) ~ 1,
#     family = binomial(link="identity"),
#     data = data)

# Check the default prior!

# Choose a meaningful prior for the Intercept parameter, use lb=0, ub=1
# Fit the model & verify convergence.
# Re-run the analysis for different priors.

# Test the hypothesis that mean survival is > 2/3.

