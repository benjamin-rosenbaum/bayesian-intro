# Practical 1: Statistical modeling and likelihood
# by Benjamin Rosenbaum

# We learn about the likelihood function and do some didactic examples of maximum 
# likelihood estimation (MLE) with the `optim()` function.

rm(list=ls())
library("sfsmisc") # mathematical integration through data points
try(dev.off())

# example 1: survival rate -----------------------------------------------------

# Data from 3 habitats and their deer population. 
# We generate population size before and after the winter.

data = data.frame(total = c(18, 25, 20),
                  survived = c(9, 10, 8))

# Q: What is the average survival rate?
# Deterministic part: theta (just fitting a mean, or the intercept if you will)
# Stochastic part:    survived ~ Binomial(total,theta)

# Naively, we could just calculate rate for each habitat and compute their mean.
# But we can do better. Use the statistical model.

data$survived/data$total
mean(data$survived/data$total)

# Probability density function p(y|theta)

# If we know the survival rate, we can compute probability for each outcome

x = 0:20
y = dbinom(x=x, size=20, prob=0.3)
barplot(y~x, xlab="Survived out of 20", ylab="Probability")
sum(y)

# same for a different value of survival rate

x = 0:20
y = dbinom(x=x, size=20, prob=0.75)
barplot(y~x, xlab="Survived out of 20", ylab="Probability")
sum(y)

# Likelihood L(theta)=p(theta|y)

# Likelihood is probability in reverse. For a given datapoint, likelihood is the
# probability of that observation as a function of parameter value 
# (survival prob. theta)

theta = seq(0,1, length.out=100)
theta

# likelihood function for parameter theta and 1st datapoint

lik = dbinom(x=9, size=18, prob=theta)
lik
plot(theta, lik, type="l", xlab="theta", ylab="Likelihood")

# likelihood function for parameter theta and 2nd datapoint

lik = dbinom(x=10, size=25, prob=theta)
lik
plot(theta, lik, type="l", xlab="theta", ylab="Likelihood")

# ATTN: Likelihood L(theta) is not a probability density function for theta.
# It does not integrate to 1. This means we can use it for maximum likelihood,
# but not for statistical inference without knowing the full integral.

integral = integrate.xy(theta, lik) 
integral

# However, if we scale (divide) L(theta) by its integral, this new function 
# integrates to 1. But this is unfeasible in most applications (>1 parameter). 

plot(theta, lik/integral, type="l", xlab="theta")
integrate.xy(theta, lik/integral) 

# OK, so far we only calculated likelihood of single datapoints. 
# Likelihood function for the whole dataset is the product of these likelihoods.

lik = dbinom(x=data$survived[1], size=data$total[1], prob=theta)*
  dbinom(x=data$survived[2], size=data$total[2], prob=theta)*
  dbinom(x=data$survived[3], size=data$total[3], prob=theta)

plot(theta, lik, type="l", xlab="theta", ylab="Likelihood")

# Maximum likelihood means finding the parameter value for which the observed 
# data is most likely to have occurred. Here we use GLM to find that value

model1 = glm( cbind(survived, total-survived) ~ 1, data=data, 
              family=binomial)
summary(model1)

# What is happening? GLM estimates a negative survival probability?
# No, the binomial family uses a log-link as default. We can override the 
# default by specifying an identity link function (parameter is estimted on its
# original scale). More on that in Lesson 5.

model1 = glm( cbind(survived, total-survived) ~ 1, data=data, 
              family=binomial(link="identity"))
summary(model1)
confint(model1)

# Here, we get a mean survival prob of 0.43 with 95%-CI [0.31, 0.55]. 
# CIs are not computed from the full likelihood, but rely on a local approximation.



# example 2: mean body size of a mammal population -----------------------------

# Now, bodysize is a continuous response. Generate the weight of 7 individuals.

data = data.frame(weight = c(104, 120, 118, 115, 99, 110, 102))
boxplot(data$weight, ylab="Weight")

# Q: What's the mean bodysize?
# Deterministic part: mu (just fitting a mean, or the intercept if you will)
# Stochastic part:    weight ~ Normal(mu,sigma)

# first with known sigma. just 1 parameter mu 
sigma = 10

# probability for given mu and sigma

x = seq(50, 150, length.out=100)
y = dnorm(x=x, mean=100, sd=sigma)
plot(x,y, type="l")

# another mean

y = dnorm(x=x, mean=110, sd=sigma)
plot(x,y, type="l")

# likelihood = probability density of data for a parameter

mu = seq(60, 130, length.out=100)

# likelihood function for 1 datapoint

lik = dnorm(data$weight[1], mean=mu, sd=sigma)
plot(mu, lik, type="l")

# likelihood function for all datapoints

lik = dnorm(data$weight[1], mean=mu, sd=sigma)*
  dnorm(data$weight[2], mean=mu, sd=sigma)*
  dnorm(data$weight[3], mean=mu, sd=sigma)*
  dnorm(data$weight[4], mean=mu, sd=sigma)*
  dnorm(data$weight[5], mean=mu, sd=sigma)*
  dnorm(data$weight[6], mean=mu, sd=sigma)*
  dnorm(data$weight[7], mean=mu, sd=sigma)
plot(mu, lik, type="l")

# In our statistical model, there is a second parameter sigma
# Likelihood is a function of both parameters L(mu,sigma)

# Write likelihood as a function but use negative log-likelihood (NLL)
# This makes things easier since L can be VERY small and numerically unstable

likelihood = function(parameters, weights){
  lik = dnorm(weights, mean=parameters[1], sd=parameters[2], log=TRUE) 
  return(-sum(lik)) 
}

# try to minimize the NLL

likelihood(c(110,10), data$weight)
likelihood(c(111,10), data$weight)
likelihood(c(109,10), data$weight)
likelihood(c(108,10), data$weight)
likelihood(c(109,10), data$weight)
likelihood(c(109,11), data$weight)
likelihood(c(109,9), data$weight)
likelihood(c(109,8), data$weight)
likelihood(c(109,7), data$weight)
likelihood(c(109,8), data$weight)
likelihood(c(110,8), data$weight)
likelihood(c(108,8), data$weight)
likelihood(c(109,8), data$weight)

# iterative algorithm to search for optimum (mu,sigma) with initial guess (100,10)

ml = optim(fn = likelihood, 
           par = c(100,10), 
           weights = data$weight) # the data
ml

# lm-solution

model2 = lm(weight ~ 1, data=data)
summary(model2)

# visualize 2-dimensional likelihood

mu.plot = seq(from=80, to=140, length.out=200)
sigma.plot = seq(from=5, to=15, length.out=200)

test = expand.grid(mu=mu.plot, sigma=sigma.plot)
test$lik = NA
for(i in 1:nrow(test)){
  test$lik[i] = likelihood(c(test$mu[i], test$sigma[i]), data$weight)
}
lik.plot = matrix(test$lik, nrow=length(mu.plot), ncol=length(sigma.plot))

image(mu.plot, sigma.plot, exp(-lik.plot),
      xlab="mu", ylab="sigma")# , col = hcl.colors(10, "terrain"))

# plot  ML solution

points(ml$par[1], ml$par[2], col="white", lwd=2)



# example 3: body size vs age --------------------------------------------------

# Now we have a continuous predictor age. Generate data for 7 individuals.

data = data.frame(weight = c(104, 120, 118, 115, 99, 110, 102),
                  age    = c(10, 12, 11, 11, 9, 11, 10))
plot(data$age, data$weight)

# Q: What's the average growth per year? (Slope in age)
# Deterministic part: mu = a+b*age
# Stochastic part:    weight ~ Normal(mu,sigma)

# We could easily fit this model with LM

model3 = lm(weight ~ age, data=data)
summary(model3)

coef(model3)
intercept = coef(model3)[1]
slope = coef(model3)[2]

# But for didactic purposes, we're doing it the hard way. Maximum likelihood.
# We code the likelihood function (probability density of all datapoints for a 
# given set of parameters a,b,sigma)

likelihood=function(parameters, weights, ages){
  lik = dnorm(weights, 
              mean=parameters[1] + parameters[2]*ages, 
              sd=rep(parameters[3],length(weights)), 
              log=TRUE) 
  return(-sum(lik)) # negative log likelihood NLL
}

# This is the LM solution and its likelihood

plot(data$age, data$weight)
abline(intercept, slope)
likelihood(c(intercept, slope, 3.25), data$weight, data$age)

# Other solutions fit the data worse and have a higher NLL
# = lower likelihood as the LM solution

abline(26, 7.5, col="grey")
likelihood(c(26, 7.5, 3.25), data$weight, data$age)

abline(35, 7.5, col="grey")
likelihood(c(35, 7.5, 3.25), data$weight, data$age)

abline(45, 6.0, col="grey")
likelihood(c(45, 6.0, 3.25), data$weight, data$age)

# iterative algorithm to search for optimum (a,b,sigma) with initial guess (100,5,8)

ml = optim(fn = likelihood, 
           par = c(100,5,8), 
           weights = data$weight, # data
           ages = data$age)       # data
ml

# That's VERY close to the LM solution

summary(model3)

# visualize 2-dimensional likelihood (intercept & slope), sigma fixed

int.plot = seq(from=15, to=35, length.out=200)
slope.plot = seq(from=5, to=10, length.out=200)

test = expand.grid(int=int.plot, slope=slope.plot)
test$lik = NA
str(test)
for(i in 1:nrow(test)){
  test$lik[i] = likelihood(c(test$int[i], test$slope[i], 5.0), data$weight, data$age)
}
lik.plot = matrix(test$lik, nrow=length(int.plot), ncol=length(slope.plot))

image(int.plot, slope.plot, exp(-lik.plot),
      xlab="Intercept", ylab="Slope")# , col = hcl.colors(10, "terrain"))

# plot  ML solution

points(ml$par[1], ml$par[2], col="white", lwd=2)

# There is a strong negative correlation between Intercept and Slope.
# This is common if the predictor is not centered (mean=0).
# In extreme cases, that can cause numerical problems and wrong ML solutions.



# example 4: body size vs age (centered) ---------------------------------------

# We repeat the same analysis, but this time the predictor is centered

data = data.frame(weight = c(104, 120, 118, 115, 99, 110, 102),
                  age    = c(10, 12, 11, 11, 9, 11, 10))
data$agecenter = data$age-mean(data$age)
plot(data$agecenter, data$weight)

model4 = lm(weight ~ agecenter, data=data)
summary(model4)

coef(model4)
intercept = coef(model4)[1]
slope = coef(model4)[2]

likelihood=function(parameters, weights, ages){
  lik = dnorm(weights, 
              mean=parameters[1] + parameters[2]*ages, 
              sd=rep(parameters[3],length(weights)), 
              log=TRUE) 
  return(-sum(lik)) 
}

plot(data$agecenter, data$weight)

abline(intercept, slope)
likelihood(c(intercept, slope, 3.25), data$weight, data$agecenter)

abline(105, 7.5, col="grey")
likelihood(c(105, 7.5, 3.25), data$weight, data$agecenter)

abline(115, 7.5, col="grey")
likelihood(c(115, 7.5, 3.25), data$weight, data$agecenter)

abline(110, 6.0, col="grey")
likelihood(c(110, 6.0, 3.25), data$weight, data$agecenter)

ml = optim(fn = likelihood, 
           par = c(100,5,8), 
           weights = data$weight, # data 
           ages = data$agecenter) # data
ml
summary(model4)

# visualize 2-dimensional likelihood (intercept & slope), sigma fixed
int.plot = seq(from=100, to=120, length.out=200)
slope.plot = seq(from=5, to=10, length.out=200)

test = expand.grid(int=int.plot, slope=slope.plot)
test$lik = NA
str(test)
for(i in 1:nrow(test)){
  test$lik[i] = likelihood(c(test$int[i], test$slope[i], 5.0), data$weight, data$agecenter)
}
lik.plot = matrix(test$lik, nrow=length(int.plot), ncol=length(slope.plot))

image(int.plot, slope.plot, exp(-lik.plot),
      xlab="Intercept", ylab="Slope")# , col = hcl.colors(10, "terrain"))

# plot  ML solution
points(ml$par[1], ml$par[2], col="white", lwd=2)

# By using a centered predictor, the correlation between intercept & slope is resolved.



# exercise: two predictors -----------------------------------------------------

# Now we have more data and a second predictor for avg annual temperature.
# Write a model that includes both predictors for LM and ML.
# Test different initial guesses for the optim function and see if they all 
# converge to the same result.

data = data.frame(weight = c(65,80,129,41,77,133,75,75,87,109,136,134,91,49,127,120,108,126,81,112),
                  age    = c(7,7,14,6,10,15,9,8,10,13,14,15,9,7,15,13,13,13,7,12),
                  temperature = c(5.9,6.9,1.5,9.6,9.0,6.9,8.0,10.2,4.8,7.6,6.2,11.2,7.3,11.4,4.1,5.1,3.7,9.5,6.4,8.3))

model5 = lm(weight ~ age+temperature, data=data)
summary(model5)

