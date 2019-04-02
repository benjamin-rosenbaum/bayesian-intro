rm(list=ls())
library(rstan)
library(coda)
library(BayesianTools)
# setwd("~/Desktop/teaching Bayes")

rstan_options(auto_write = TRUE)
options(mc.cores = 3) # number of CPU cores

# in this session, we will learn how to fit a model and to interpret the output. 
# specifically, we learn how to deal with the posterior distribution to make inference and predictions.
# Again, we will use linear regression.

# statistical model: linear regression
# y_i ~ normal(\mu_i, \sigma)
# mu_i = a+b*x_i

#------------------------------------------------------------------------------
# generate data 
#------------------------------------------------------------------------------

set.seed(123) # initiate random number generator for reproducability

n=50

a=1.0
b=0.5
c=0.4
sigma=0.15

x = runif(n=n, min=0, max=1)
y = rnorm(n=n, mean=a+b*x+c*x^2, sd=sigma)

df = data.frame(x=x,
                y=y)

plot(df)


#------------------------------------------------------------------------------
# define stan model
#------------------------------------------------------------------------------
stan_code = '
data {
  int n;
  vector[n] x;
  vector[n] y;
}
parameters {
  real a;  
  real b;  
  real<lower=0> sigma;  // standard deviation
}
model {
  // priors
  a ~ normal(0, 10);
  b ~ normal(0, 10);
  sigma ~ normal(0, 10);
  // likelihood
  y ~ normal(a+b*x, sigma);
}
'

#------------------------------------------------------------------------------
# prepare data and compile model
#------------------------------------------------------------------------------
data = list(n=n, 
            x=df$x, 
            y=df$y)

stan_model = stan_model(model_code=stan_code)
# save(file="stan_model_test.RData", list="stan_model")
# load("stan_model_test.RData")

#------------------------------------------------------------------------------
# MCMC sampling 
#------------------------------------------------------------------------------
fit  = sampling(stan_model,
                data=data,
                chains=3,
                iter=2000,
                warmup=1000
)

#------------------------------------------------------------------------------
# explore the posterior distribution
#------------------------------------------------------------------------------

# first, we look at the output and check n_eff and Rhat.
# Rhat<1.01, that's good. n_eff looks good, too (compare to n_total)

print(fit, digits=3, probs=c(0.025, 0.975))

# plot the standard Stan output. 

plot(fit)

# plotting from coda package (covert fit object to mcmc.list object)
# also shows traceplots of the 3 chains.
# They look like a fat hairy caterpillar, 
# so we assume the chains are a good representation of the true posterior distribution.

plot(As.mcmc.list(fit)) # from coda package

# What is the posterior exactly? 
# Convert fit object into a matrix. 
# Now the 3 chains are concatenated to 1.
# It has n_total=3000 rows and 1 column per parameter 
# (3 parameters + "lp__", we ignore the last one)

posterior=as.matrix(fit)
str(posterior)
head(posterior)

# each column contains all posterior samples of 1 parameter.
# We can look at this (marginal) posterior distribution. 
# This is the same as the plotting commands above.
# We can index the columns by number or by name.

str(posterior[, 1])
head(posterior[, 1])
hist(posterior[, 1])
hist(posterior[, "a"])

# each row contains one sample of the multidimensional posterior.
# important: each draw / sample consists of a multidimensional vector for a,b,sigma.
# if you change the order of one column (permutation), 
# the whole thing is not a representation of the posterior anymore!
# each element of the matrix is linked to the other elements in that row!

str(posterior[1, ])
posterior[1, ]

# the reason is that the parameters can be correlated. 
# typically, there is some correlation between intercept and slope in linear regression.
# (especially if the data is not centered).
# correlationPlot() shows pairwise plots and correlation coefficients.
# some correlation is generally not a problem in Bayesian statistics! 
# Perfect correlation (samples perfectly distributed along a line or a curve), 
# however, would indicate some problem with the model (unidentifiability).

pairs(fit, pars=c("a","b","sigma"))
correlationPlot(posterior[, 1:3], thin=1) # from BayesianTools package

# Now it's time for some inference!
# We want to test if there is a positive effect of predictor x on response y.
# Posterior samples of slope b represent posterior distribution of b given the data p(b|y).
# So we can actually compute the posterior probability of a positive effect given the data P(b>0|y)
# How to do that? just count the number the event (b>0) occurs in the posterior samples. 
# divide by total number of samples and that's the probability!
# Given the data, we are 100% sure that the effect is positive.

hist(posterior[, "b"], xlim=c(0,max(posterior[, "b"])))
abline(v=0, col="red", lwd=2)
sum(posterior[, "b"]>0)/nrow(posterior)

# Similarly, we can compute the probability of the effect being larger than 1
# or effect being in the interval [0.9, 1.1]

sum(posterior[, "b"]>1)/nrow(posterior)
sum( (posterior[, "b"]>0.9) & (posterior[, "b"]<1.1) )/nrow(posterior)

#------------------------------------------------------------------------------
# posterior predictions 
#------------------------------------------------------------------------------

# Usually, it is *not* sufficient just to check if the MCMC sampler converged.
# That doesn't tell us anything about if our statistical model 
# (deterministic part and stochastic part) describes the data adequately!
# For that, we have to compare observed and predicted values

# Each row of the posterior matrix contains a sample of the posterior, 
# i.e. intercept a and slope b.
# We can evaluate or plot the deterministic model (regression line a+b*x) 
# using these parameters.
# E.g. plot the deterministic model for the first sample using the abline() command.

plot(df)
abline(posterior[1,"a"], posterior[1,"b"], col=adjustcolor("red", alpha.f=0.3))

# Or plot the deterministic model for the first 100 samples. 
# We can see the uncertainty associated with the predictions.
# Remember that each row is a sample, i.e. we have to use intercept a_i and slope b_i.
# Never mix up the order, a_i with slope b_j (i not equal j)!

plot(df)
for(i in 1:100){
  abline(posterior[i,"a"], posterior[i,"b"], col=adjustcolor("red", alpha.f=0.3))
}

# abline() is a fancy command for plotting lines, 
# but if we want a more generalized approach (more complex models later), 
# we have to code the deterministic model ourself.
# x.pred is a vector of predictor values for which we want to make predictions
# y.cred is a matrix that will contain all predictions. 
# We call it "cred" because we will use it for computing 
# "credible intervals" / "confidence intervals" of the deterministic model.
# See below for "prediction intervals"
# In Bayesian statistics, everything is a distribution. 
# So also the predictions will be a distribution. 
# There are 3000 samples in the posterior, 
# i.e. 3000 parameter combinations of intercept and slope.
# This means we can make 3000 predictions and these will be 
# samples from a posterior predictive distribution.

x.pred = seq(from=0, to=1, by=0.1)
y.cred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))

for(i in 1:nrow(posterior)){
  y.cred[i, ] = posterior[i,"a"] + posterior[i,"b"]*x.pred
}

# The element y.cred[i,j] (ith row, jth column) is the prediction for 
# sample i (using parameters a_i,b_i) for predictor value x_j
# Each row i contains predictions for all x using single parameter set a_i,b_i.
# Each column j contains 3000 predictions (for all posterior samples) for one predictor value x_j.

head(y.cred)
hist(y.cred[,1])

# As above, we can plot the first 100 predictions

plot(df)
for(i in 1:100){
  lines(x.pred, y.cred[i, ], col=adjustcolor("red", alpha.f=0.3))
}

# Since each column contains samples of a posterior distribution, 
# we can make statistics, e.g. mean or confidence intervals. 
# In Bayesian stats, these confidence intervals are often called "credible intervals". 
# We will now plot the mean and the 90% credible intervals (using 5% and 95% quantiles). 
# Why 90% and not 95%? 95% is an arbitrary number. Choose your own credible interval!
# We use the apply() function to use the mean() and quantile() commands on each column of the matrix.

plot(df)

y.cred.mean = apply(y.cred, 2, function(x) mean(x)) 
lines(x.pred, y.cred.mean, col="red", lwd=2)

y.cred.q05 = apply(y.cred, 2, function(x) quantile(x, probs=0.05)) 
lines(x.pred, y.cred.q05, col="red", lwd=2, lty=2)

y.cred.q95 = apply(y.cred, 2, function(x) quantile(x, probs=0.95)) 
lines(x.pred, y.cred.q95, col="red", lwd=2, lty=2)

# A statistical model contains a deterministic and stochastic part. 
# The credible / confidence intervals are computed using the 
# distribution of the deterministic part only!
# They are confidence intervals for the regression line, not for the data!
# Now we will compute true prediction intervals also using the stochastic model part 
# (data are normally distributed around regression line with standard deviation \sigma).

# y.pred is structured as y.cred above:
# Each row i contains predictions for all x using 
# a single parameter set a_i,b_i (and sigma_i).
# Each column j contains 3000 predictions (for all posterior samples)
# for one predictor value x_j.
# But now, each prediction is a random draw from normal(a_i+b_i*x,sigma_i)
# (deterministic part a_i+b_i*x was already computed in y.cred)

y.pred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))

for(i in 1:nrow(posterior)){
  y.pred[i, ] = rnorm(n=length(x.pred), mean=y.cred[i, ], sd=rep(posterior[i, "sigma"],length(x.pred)) )
}

plot(df)

lines(x.pred, y.cred.mean, col="red", lwd=2)
lines(x.pred, y.cred.q05, col="red", lwd=2, lty=2)
lines(x.pred, y.cred.q95, col="red", lwd=2, lty=2)

y.pred.mean = apply(y.pred, 2, function(x) mean(x)) 
lines(x.pred, y.pred.mean, col="blue", lwd=2)

y.pred.q05 = apply(y.pred, 2, function(x) quantile(x, probs=0.05)) 
lines(x.pred, y.pred.q05, col="blue", lwd=2, lty=2)

y.pred.q95 = apply(y.pred, 2, function(x) quantile(x, probs=0.95)) 
lines(x.pred, y.pred.q95, col="blue", lwd=2, lty=2)

legend("topleft", legend=c("90% credible","90% prediction"), lwd=c(2,2), col=c("red","blue"), bty="n", lty=c(2,2))

# The 90% credible interval (red) tells us that we are 
# 90% sure that the regression line is in that interval.
# The 90% prediction interval (blue) tells us that we are 
# 90% sure that the data are in that interval.

# 4 out of 50 datapoints are outside the prediction interval. 
# 46 out of 50 datapoints are inside the prediction interval, that's 92%.

# Note: you can also make predictions while fitting using the generated_quantities{} block.

#------------------------------------------------------------------------------
# posterior predictions - obs vs pred for data
#------------------------------------------------------------------------------

# In the section above we used a sequence of predictor values x.pred for making nice plots.
# For model validation, we should make predictions for the actual data.
# Then we can compare observed and predicted values 
# (even if we have many predictors / groups and nice plots aren't possible).

x.pred = df$x
y.cred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))

for(i in 1:nrow(posterior)){
  y.cred[i, ] = posterior[i,"a"] + posterior[i,"b"]*x.pred
}

y.cred.mean = apply(y.cred, 2, function(x) mean(x))
y.cred.q05 = apply(y.cred, 2, function(x) quantile(x, probs=0.05))
y.cred.q95 = apply(y.cred, 2, function(x) quantile(x, probs=0.95))

y.pred = matrix(0, nrow=nrow(posterior), ncol=length(x.pred))
for(i in 1:nrow(posterior)){
  y.pred[i, ] = rnorm(n=length(x.pred), mean=y.cred[i, ], sd=rep(posterior[i, "sigma"],length(x.pred)) )
}

y.pred.mean = apply(y.pred, 2, function(x) mean(x)) 
y.pred.q05 = apply(y.pred, 2, function(x) quantile(x, probs=0.05)) 
y.pred.q95 = apply(y.pred, 2, function(x) quantile(x, probs=0.95)) 

plot(df$y, y.pred.mean, ylim=c(min(df$y), max(df$y)), xlab="observed", ylab="predicted")
abline(0,1)
for (i in 1:n){
  lines(c(df$y[i], df$y[i]), c(y.pred.q05[i], y.pred.q95[i]))
}

# residual plots

plot(df$y, df$y-y.pred.mean)
abline(0,0)

# There seems to be a little systematic underfitting going on for high observed values. 
# (the dataset contains a small quadratic effect, see above)

# Note: you can also make predictions while fitting using the generated_quantities{} block.

#------------------------------------------------------------------------------
# model comparison, AIC
#------------------------------------------------------------------------------
# will not be covered here. 
# use the "loo" package for an information criterion that you can use similar to AIC.
# you have to calculate the pointwise log-likelihood values in your model in the "generated quantities{}" block.
# see https://cran.r-project.org/web/packages/loo/vignettes/loo2-with-rstan.html

#------------------------------------------------------------------------------
# pitfall of predictions
#------------------------------------------------------------------------------

# That's a lot of code above just for predictions. 
# Can't we just use the mean fitted parameters to make predictions?
# Please **never** do that! In Bayesian stats, everything is a distribution, also the predictions. 
# For simple linear models, the mean of predictions can be equal to the predictions using the mean parameters.
# This is not the case for more complex models!

print(fit)
summary.fit = summary(fit)$summary

summary.fit

plot(df)
abline(summary.fit["a","mean"], summary.fit["b","mean"], col="red", lwd=2)
points(x.pred,y.cred.mean, col="red", pch="+")

legend("topleft", legend=c("prediction using mean parameters (WRONG!)","means of predictions"), bty="n",
       lty=c(1,NA), pch=c(NA,"+"), col=c("red","red"), lwd=c(2,2))

