rm(list=ls())
library(rstan)
library(coda)

set.seed(123) # initiate random number generator for reproducability
 
#------------------------------------------------------------------------------
# generate data
#------------------------------------------------------------------------------
n=100

a=1
b=2
sigma=0.5

x = runif(n=n, min=0, max=1)
y = rnorm(n=n, mean=a+b*x, sd=sigma)

df = data.frame(x=x,
                y=y)

plot(df)

#------------------------------------------------------------------------------
# 2 statistical model, deterministic and stochastic part
#------------------------------------------------------------------------------

# statistical model: linear regression
# y_i ~ normal(\mu_i, \sigma)
# mu_i = a+b*x_i

#------------------------------------------------------------------------------
# define stan model
#------------------------------------------------------------------------------

# The Stan code can be written in the R script and is saved as text using '...'.
# We use the 3 blocks data{}, parameters{} and model{}. 
# The first 2 blocks just contain the definition of all variables. 
# Note that sigma is constrained to be positive. 
# In the model block we define priors for all parameters and the likelihood function, 
# which contains a deterministic part (a+b*x) 
# and a stochastic part (y~normal(..., sigma).)
# We can use '//' to comment the code (same as '#' in R).
# The likelihood definition   y ~ normal(a+b*x, sigma); is short for
# for (i in 1:n){
#   y[i] ~ normal(a+b*x[i], sigma);
# }
# We use (very) weakly informative priors 
# (normal distribution with zero mean and sd of 10). 
# Notice that sigma has a lower boundary of 0, 
# i.e. the prior for sigma is a positive half-normal distribution.
# We use a general expectation as prior information (e.g. intercept of 100 wouldn't make sense), 
# but not summary statistics of the data (e.g. slope from lm()-fit etc. - please never do that!)
# More on priors tomorrow!

stan_code = '
data {
  int n;
  vector[n] x;
  vector[n] y;
}
parameters {
  real a; // intercept  
  real b; // slope
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
# prepare data
#------------------------------------------------------------------------------

# The data have to be prepared in R according to the data{} block in the Stan code. 
# It is saved as a named list.

data = list(n=n, 
            x=df$x, 
            y=df$y)
data

#------------------------------------------------------------------------------
# compile model (translates to fast C++ code)
#------------------------------------------------------------------------------
# usually, we run more than 1 MCMC chain. 
# These calculations are totally independent from each other and can be run simultaneously. 
# detectCores() shows the number of available CPU cores. 

parallel::detectCores(all.tests = FALSE, logical = TRUE)

# Now we tell Stan how many CPU cores we want to use.
# We want to run 3 MCMC chains in parallel, so we set this number to 3.

rstan_options(auto_write = TRUE)
options(mc.cores = 3) # number of CPU cores used

# Until yet, the Stan code is just saved as text.
# stan_model() compiles this text into fast C++ code and we save this as an object. 
# This can take a minute!

stan_model = stan_model(model_code=stan_code)
# save(stan_model, file="stan_model.RData")
# load("stan_model.RData")

#------------------------------------------------------------------------------
# MCMC sampling 
#------------------------------------------------------------------------------

# Now we can start the MCMC sampler using sampling(). 
# We hand over the compiled object and the data. 
# Additionally, we tell the sampler how many Markov chains to run, 
# how many iterations per chain, and how many of these are just used for warmup 
# (discarded afterwards).
# The total number of samples will be chains*(iter-warmup).
# More options see help function ?sampling. 
# One important option we didn't use is specifying initial values for the chains. 
# If it is unspecified, Stan uses random initial values from the interval [-2,2].

fit  = sampling(stan_model,
                data=data,
                chains=3,
                iter=2000,
                warmup=1000
)

# Printing the results gives is information on the marginal (i.e. 1d / single parameter) estimates 
# (remember, the posterior distribution is multidimensional):
# mean, standard deviation and quantiles / confidence intervals. 
# Convergence of the MCMC sampling is monitored by n_eff and Rhat. 
# n_eff is the number of effective (uncorrelated) samples of the posterior 
# and is generally lower than n_total. 
# But if n_eff was only a small fraction of n_total, 
# i.e. 1%, this would indicate some problems with our model.
# Rhat checks if all chains behave similarly, values close to one (<1.1) indicate good behavior.
# The last "parameter" is the log posterior denstity, we will ignore it for now.

print(fit)

# You can also specify the output, printing only some parameters 
# or only 95% confidence intervals. 
# Btw, the choice of 95% confidence interval is rather arbitrary. 
# We can also print a 90% confidence interval (5% cut off at both sides). 

print(fit,
      digits=3,
      probs=c(0.05, 0.95),
      pars=c("a","b","sigma")
      )

# The posterior distribution can also be plotted.

plot(fit)

# Or just for certain parameters:

plot(fit, pars=c("a","b"))

# There is another option from the coda package. 
# First we transform the stan object

posterior = As.mcmc.list(fit)

# Plotting this object shows us traceplots of the chains 
# and marginal density functions of the posterior distribution.
# The chains are plotted in different colours. 
# We also can visually inspect if all chains behave similarly.  
# They should look like a "fat hairy caterpillar".

plot(posterior)

# You can also specify which parameters to plot, either indicing by name or number.

plot(posterior[, c("a", "b")])
plot(posterior[, 1:2])

# The posterior distribution is multidimensional (here, d=2). 
# So visualization of the whole posterior is limited, 
# but we can at least look at the pairwise 2d marginals

pairs(fit, pars=c("a","b","sigma"))

# nicer pairs plot from BayesianTools package, 
# also computes the pairwise correlation coefficients of the parameters
library(BayesianTools)
correlationPlot(as.matrix(fit)[,1:3], thin=1)

#------------------------------------------------------------------------------
# other stan plotting options 
#------------------------------------------------------------------------------
stan_plot(fit)
stan_trace(fit)
stan_hist(fit)
stan_dens(fit)
stan_scat(fit, pars=c("a","b"))


#------------------------------------------------------------------------------
# just for demonstration: maximum a posteriori estimation 
#------------------------------------------------------------------------------

# MCMC sampling globally explores the posterior distribution. 
# Stan can also find the paramter combination that maximizes the posterior probability density function.
# Similiar to maximum likelihood estimation (MLE), this is called maximum a posteriori (MAP).
# It just computes a point estimate, without confidence intervals / uncertainty / sd.

map = optimizing(stan_model,
                 data=data
)

print(map)

