rm(list=ls())
library(manipulate)

set.seed(123) # initiate random number generator for reproducability

#------------------------------------------------------------------------------
# 1 generate data
#------------------------------------------------------------------------------
n=50

y = rnorm(n=n, mean=1.0, sd=2.0)

# examine data:
y
hist(y)
points(y, jitter(rep(0, times=n), factor=10))

#------------------------------------------------------------------------------
# 2 statistical model, deterministic and stochastic part
#------------------------------------------------------------------------------

# statistical model: estimate mean and standard deviation
# y_i ~ normal(\mu, \sigma)

#------------------------------------------------------------------------------
# 3 the likelihood function 
#------------------------------------------------------------------------------

# likelihood function for a single datapoint probability density function p(y_i|mu,sigma)= ... exp(...)
# single data point, likelihood function L can be computed with dnorm() function (for a given parameter combination mu, sigma)

# likelihood of first data point:
i = 1
L = dnorm(x=y[i], mean=0, sd=1)
y[i]
L

# plot
curve(dnorm(x, mean=0, sd=1), from=-4, to=6, xlab="y", ylab="p(y)", ylim=c(0, 0.6))
points(y[i],0)
lines(c(y[i],y[i]), c(0, L), col="red")
lines(c(-50,y[i]), c(L, L), col="red")

# all datapoints

# dnorm can calculate all p(y_i|mu,sigma) for given parameter combination mu,sigma at once
L = dnorm(x=y, mean=0, sd=1)

# L is a vector now
L

# the likelihood function of all datapoints for a given parameter combination mu,sigma is the product of all single values
# p(y_1,...y_n|mu,sigma)=p(y_1|mu,sigma)*...*p(y_n|mu,sigma)
# This holds because observations are independent
prod(L)

# There is a problem. each p(y_i|mu,sigma)<1. So multiplying them all results in an extremely small number. 


# Better: use log p(y_1,...y_n|mu,sigma)  = log( p(y_1|mu,sigma)*...*p(y_n|mu,sigma)  )) = log p(y_1|mu,sigma) + ... + log p(y_n|mu,sigma)
# log of a product is equal to sum of logs. We minimize the negative log likelihood (NLL) to find model parameters, which is equivalent to maximum likelihood (mathematical convention is minimization instead of maximization)
-sum(log(L))

# visualize likelihood function of all datapoints for given parameters \mu and \sigma
# you can play around with mu and sigma. which combination maximizes likelihood of all datapoints at once?
curve.data <- function(mean, sd)
{
  curve(dnorm(x, mean=mean, sd=sd), from=-4, to=6, xlab="y", ylab="p(y)", ylim=c(0, 0.6))
  for(i in 1:n){
    lines(c(y[i],y[i]),c(0,dnorm(x=y[i], mean=mean, sd=sd)))
  }
}

manipulate(curve.data(mean, sd), 
           mean=slider(-3, 3, step=0.1, initial=0), 
           sd=slider(0.1,4, step=0.1, initial=1) )


#------------------------------------------------------------------------------
# 4 maximum likelihood with optim()
#------------------------------------------------------------------------------

# We can use mathematical algorithms to search for the best parameter combination automatically
# First, we define a function that directly calculates the NLL for the data and a given parameter combination

nll.function = function(data, par){ # nll: negative log likelihood
  LL = dnorm(x=data, mean=par[1], sd=par[2], log=TRUE) # LL: log-likelihood
  NLL = -sum(LL)
  return(NLL)
}

# example: for mu=0, sd=1

nll.function(data=y, par=c(0.0, 1.0))

# optim() function automatically searches for parameters that minimize the NLL

optim(par=c(0.0, 1.0), # initial guess mu=0, sd=1
      fn=nll.function,
      data=y)

