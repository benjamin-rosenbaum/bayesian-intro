# Practical 7: Introduction to Stan
# by Benjamin Rosenbaum

# We code some classical models in Stan, using the rstan package. 
# We learn how to interpret and assess Stan model output. 
# Posterior predictions unfortunately have to be coded manually.

rm(list=ls())
library("rstan")
library("loo")
library("ecostats")
library("Data4Ecologists")
library("bayesplot")
library("ggplot2")
library("arm")
library("RColorBrewer")

setwd("~/Nextcloud/Teaching brms/Practical_07")

# linear regression ------------------------------------------------------------

# global plants dataset, log(height) ~ latitude
# scale predictor latitude

data("globalPlants")
data = globalPlants
plot(data$lat, log(data$height), xlab="Latitude", ylab="log(height)")

# Put data in named list for Stan. Names must be identical to Stan code names

stan.data = list(N = nrow(data),
                 x = as.vector(scale(data$lat)), # convert to vector to get rid of additional infos from scale function
                 y = log(data$height))

# Compiling code & MCMC with single function stan()

fit1 = stan(file="lm.stan", data=stan.data)

# check convergence 

print(fit1)
print(fit1, probs=c(0.05, 0.95))

# can extract summary statistics in dataframe / matrix if needed

summ.table = summary(fit1, probs=c(0.05, 0.95))$summary 
print(summ.table, digits=3)

# from bayesplot package. 
# alternatively, can also use stan_trace, stan_hist, stan_dens ... from rstan 

mcmc_hist(fit1)
mcmc_trace(fit1)

# side-by-side as in brms

mcmc_combo(fit1)
mcmc_combo(fit1, combo=c("hist","trace"), pars=c("a","b","sigma"))

# posterior predictions (~ conditional_effects)

# extract posterior samples as matrix or dataframe.
# each column = 1 parameter, each row = 1 posterior sample

post = as.matrix(fit1)
head(post[, 1:3])

# set range of predictor values for plotting. here only 1 predictor. 
# for multiple predictors, choose fixed level for other predictor(s)

xmin = min(stan.data$x)
xmax = max(stan.data$x)
x.pred = seq(xmin, xmax, length.out=100)

# fitted (deterministic part)

# columns of this matrix correspond to x.pred values,
# each row contains a regression line

y.fit = matrix(NA, nrow=nrow(post), ncol=length(x.pred) )
for(i in 1:nrow(post)){
  y.fit[i, ] = post[i,"a"] + post[i,"b"]*x.pred
}
str(y.fit)

# spaghetti plot

plot(stan.data$x, stan.data$y, xlab="lat.z", ylab="log(height)", cex=0.8)
draws = sample(1:nrow(post), size=50)
for(i in draws){
  lines(x.pred, y.fit[i, ], col=adjustcolor("blue", alpha.f=0.33))
}

# extract mean and CI. summarize each column: apply(..., margin=2, function() ... )
# for fitted plot

y.fit.mean = apply(y.fit, 2, function(x) mean(x)) 
y.fit.q05 = apply(y.fit, 2, function(x) quantile(x, probs=0.05)) 
y.fit.q95 = apply(y.fit, 2, function(x) quantile(x, probs=0.95)) 

plot(stan.data$x, stan.data$y, xlab="lat.z", ylab="log(height)", cex=0.8)
polygon(c(x.pred, rev(x.pred)),
        c(y.fit.q05, rev(y.fit.q95)),
        border = NA,
        col = adjustcolor("blue", alpha.f=0.25))
lines(x.pred, y.fit.mean, col="blue", lwd=2)

# predicted (deterministic & stochastic part)

# columns of this matrix correspond to x.pred values,
# we predict data based on statistical model (normal distr. around fitted value)

y.pred = matrix(NA, nrow=nrow(post), ncol=length(x.pred) )
for(i in 1:nrow(post)){
  y.pred[i, ] = rnorm(n = length(x.pred), 
                      mean = y.fit[i, ], 
                      sd = post[i, "sigma"] )
}

# extract mean and prediction intervals for plot

y.pred.mean = apply(y.pred, 2, function(x) mean(x)) 
y.pred.q05 = apply(y.pred, 2, function(x) quantile(x, probs=0.05)) 
y.pred.q95 = apply(y.pred, 2, function(x) quantile(x, probs=0.95)) 

plot(stan.data$x, stan.data$y, xlab="lat.z", ylab="log(height)", cex=0.8)
polygon(c(x.pred, rev(x.pred)),
        c(y.pred.q05, rev(y.pred.q95)),
        border = NA,
        col = adjustcolor("red", alpha.f=0.25))
lines(x.pred, y.fit.mean, col="red", lwd=2)

# posterior predictions (on data level)

# same code as before, but use actual observations for prediction

x.pred = stan.data$x
y.fit  = matrix(NA, nrow=nrow(post), ncol=length(x.pred) )
y.pred = matrix(NA, nrow=nrow(post), ncol=length(x.pred) )
for(i in 1:nrow(post)){
  y.fit[i, ]  = post[i,"a"] + post[i,"b"]*x.pred
  y.pred[i, ] = rnorm(n=length(x.pred),
                      mean = y.fit[i, ],
                      sd = post[i,"sigma"]) 
}

# Residuals plot

y.fit.mean = apply(y.fit, 2, function(x) mean(x)) 
residuals = stan.data$y - y.fit.mean
plot(y.fit.mean, residuals, xlab="Fitted", ylab="Residuals")
abline(0,0)

# Posterior predictive checks: 
# histograms / densities of 100 predicted datasets vs true dataset

dens.data = density(stan.data$y)
plot(dens.data$x, dens.data$y, type="n", xlab="y", ylab="density")
draws = sample(1:nrow(post), size=100)
for(i in draws){
  dens = density(y.pred[i, ])
  lines(dens$x, dens$y, col=adjustcolor("grey", alpha.f=0.33))
}
lines(dens.data$x, dens.data$y, col="blue", lwd=2)

# GLM --------------------------------------------------------------------------

# Binomial regression: proportion of disinfected mice vs parasite dose
# deterministic: logit(mu)=a+b*dose  -->  mu=invlogit(a+b*dose)
# stochastic:    y ~ Binomial(total,mu)

df = read.csv("https://raw.githubusercontent.com/songsqian/eesR/refs/heads/master/R/Data/cryptoDATA.csv",
              sep=" ")
data = subset(df, Source=="Finch")

plot(data$Dose, data$Y/data$N, xlab="Dose", ylab="Proportion")

stan.data = list(N = nrow(data),
                 dose = as.vector(scale(data$Dose)), # get rid of additional infos from scale function
                 total = as.integer(data$N),
                 y = as.integer(data$Y))

fit2 = stan(file="glm.stan", data=stan.data)

print(fit2, probs=c(0.05, 0.95))
mcmc_combo(fit2, combo=c("hist","trace"), pars=c("a","b"))

# posterior predictions (~ conditional_effects)

post = as.matrix(fit2)
head(post[, 1:2])

xmin = min(stan.data$dose)
xmax = max(stan.data$dose)
x.pred = seq(xmin, xmax, length.out=100)

# fitted (deterministic part)

y.fit = matrix(NA, nrow=nrow(post), ncol=length(x.pred) )
for(i in 1:nrow(post)){
  y.fit[i, ] = invlogit(post[i,"a"] + post[i,"b"]*x.pred)
}
str(y.fit)

# spaghetti plot

plot(stan.data$dose, stan.data$y/stan.data$total, xlab="Dose", ylab="Infected", cex=0.8)
draws = sample(1:nrow(post), size=50)
for(i in draws){
  lines(x.pred, y.fit[i, ], col=adjustcolor("blue", alpha.f=0.33))
}

# extract mean and CI  for fitted plot

y.fit.mean = apply(y.fit, 2, function(x) mean(x)) 
y.fit.q05 = apply(y.fit, 2, function(x) quantile(x, probs=0.05)) 
y.fit.q95 = apply(y.fit, 2, function(x) quantile(x, probs=0.95)) 

plot(stan.data$dose, stan.data$y/stan.data$total, xlab="Dose", ylab="Infected", cex=0.8, main="Fitted")
polygon(c(x.pred, rev(x.pred)),
        c(y.fit.q05, rev(y.fit.q95)),
        border = NA,
        col = adjustcolor("blue", alpha.f=0.25))
lines(x.pred, y.fit.mean, col="blue", lwd=2)

# Binned residuals plot

x.pred = stan.data$dose
y.fit = matrix(NA, nrow=nrow(post), ncol=length(x.pred) )
for(i in 1:nrow(post)){
  y.fit[i, ] = invlogit(post[i,"a"] + post[i,"b"]*x.pred)
}
str(y.fit)
y.fit.mean = apply(y.fit, 2, function(x) mean(x)) 
y.fit.response = y.fit.mean*stan.data$total # scale up from proportion to response
y.residuals = y.fit.response-stan.data$y 

plot(y.fit.mean, y.residuals)
binnedplot(y.fit.mean, y.residuals)

# ANOVA --------------------------------------------------------------------------

# categorical predictor: species richness vs. landscape type
# deterministic: mu_i = b(landscape_i)       i=1:N
# stochastic:    S_i ~ normal(mu_i, sigma)   i=1:N
# here not dummy-coded, but effect coded:
# group-means b1,b2,b3,b4

data(birds)
data = birds[, c("S", "landscape", "log.area.")] # use area later
data$landscape = as.factor(data$landscape)
data = data[complete.cases(data), ]

plot(data$landscape, data$S, xlab="Landscape", ylab="Species richness S")

# convert factor to integer for Stan
stan.data = list(N = nrow(data),
                 M = 4,
                 y = data$S,
                 group = as.integer(data$landscape))

fit3 = stan(file="lm_anova.stan", data=stan.data)

print(fit3, probs=c(0.05, 0.95))
mcmc_combo(fit3, combo=c("hist","trace"))

# extract means and CIs directly from summary for plotting fitted

summ.table = summary(fit3, probs=c(0.05, 0.95))$summary |> round(3)
summ.table

plot(jitter(stan.data$group), stan.data$y, 
     xlim=c(0.5,4.5),
     pch=16,
     col=adjustcolor(1, alpha.f=0.33),
     xlab="landscape", ylab="S")
for(i in 1:4){
  points(i, summ.table[i, "mean"], col="blue", pch=16, cex=1.5)
  lines(c(i, i), summ.table[i, 4:5], col="blue", lwd=3)
}

# LMM --------------------------------------------------------------------------

# same data as ANOVA. ANOVA uses no pooling for landscape type
# LMM uses partial pooling for landscape type
#
# want to estimate overall mean S, but acknowledge non-independence of residuals
# here this would make sense, if actual sites of observations are not mixed, 
# but 1 geographical cluster per landscape type with multiple sites (1 site = 1 obs)
# --> spatial autocorrelation
#
# deterministic: mu_i = b(landscape_i)       i=1:N
# stochastic:    S_i ~ normal(mu_i, sigma)   i=1:N
# hierarchical:  b_j ~ normal(mu_b, sigma_b) j=1:4
# so random effects are actual group means, not deviation from grand mean mu_b

fit4 = stan(file="lmm_intercepts.stan", data=stan.data)

print(fit4, probs=c(0.05, 0.95))
mcmc_combo(fit4, combo=c("hist","trace"))

# plot fitted by extracting means and CIs from the summary

summ.table = summary(fit4, probs=c(0.05, 0.95))$summary |> round(3)
summ.table

plot(jitter(stan.data$group), stan.data$y, 
     xlim=c(0,4.5),
     pch=16,
     col=adjustcolor(1, alpha.f=0.33),
     xlab="landscape", ylab="S")
for(i in 1:4){
  points(i, summ.table[i+2, "mean"], col="blue", pch=16, cex=1.5)
  lines(c(i, i), summ.table[i+2, 4:5], col="blue", lwd=3)
}
points(0, summ.table[1, "mean"], col="purple", pch=16, cex=1.5)
lines(c(0,0), summ.table[1, 4:5], col="purple", lwd=3)

# exercise: ANCOVA S~area+landscape --------------------------------------------

# add a continuous predictor area (ANCOVA), `S ~ landscape + log.area`

# BONUS: save pointwise loglik 
# refit ANOVA with generated quantities
# refit ANCOVA with generated quantities
# LOO model comparison

# https://mc-stan.org/loo/articles/loo2-with-rstan.html

data(birds)
data = birds[, c("S", "landscape", "log.area.")]
data$landscape = as.factor(data$landscape)
data = data[complete.cases(data), ]
plot(data$log.area., data$S)
