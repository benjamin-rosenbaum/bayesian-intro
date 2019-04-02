rm(list=ls())
library(rstan)
library(coda)

rstan_options(auto_write = TRUE)
options(mc.cores = 3) 

#------------------------------------------------------------------------------
# exercise 1 - model parameters as vector
#------------------------------------------------------------------------------

# defining the vector b in the parameter block:   
# vector[2] b;  
# when defining the prior, b ~ normal(0, 10); is now short version of:
# b[1] ~ normal(0, 10);
# b[2] ~ normal(0, 10);
# (can use both)

set.seed(123) # initiate random number generator for reproducability

n=100

a=1
b=2
sigma=0.5

x = runif(n=n, min=0, max=1)
y = rnorm(n=n, mean=a+b*x, sd=sigma)

df = data.frame(x=x,
                y=y)

plot(df)

stan_code = '
data {
  int n;
  vector[n] x;
  vector[n] y;
}
parameters {
  vector[2] b;  // real b[2];
  real<lower=0> sigma;  
}
model {
  // priors
  b ~ normal(0, 10);
  
  // b[1] ~ normal(0, 10);
  // b[2] ~ normal(0, 10);

  sigma ~ normal(0, 10);
  // likelihood
  y ~ normal(b[1]+b[2]*x, sigma);
}
'

data = list(n=n, 
            x=df$x, 
            y=df$y)

stan_model = stan_model(model_code=stan_code)
# save(stan_model, file="stan_model.RData")
# load("stan_model.RData")

fit  = sampling(stan_model,
                data=data,
                chains=3,
                iter=2000,
                warmup=1000
                )

print(fit, probs=c(0.025, 0.975))

posterior = As.mcmc.list(fit)
plot(posterior[, c("b[1]", "b[2]", "sigma")])

#------------------------------------------------------------------------------
# exercise 2 - quadratic stan model
#------------------------------------------------------------------------------

# define vector[3] b;  in the parameters block containing intercept, 
# effects of linear and quadratic term. 
# again, in the prior definition b ~ normal(0, 10); is short for 
# b[1] ~ normal(0, 10);
# b[2] ~ normal(0, 10);
# b[3] ~ normal(0, 10);
# (can use both)
# Attention: we can't use the formulation "b[3]*x^2" or "b[3]*x*x" 
# because these operations aren't defined for a vector x in Stan.
# Either use the pointwise vector operation x .* x or a for loop!

stan_code_2 = '
data {
  int n;
  vector[n] x;
  vector[n] y;
}
parameters {
  vector[3] b;  
  real<lower=0> sigma;  
}
model {
  // priors
  b ~ normal({0.0, 0.0 , 0.0}, {10.0, 10.0, 10.0});
  sigma ~ normal(0, 10);
  // likelihood
  // for(i in 1:n){
  //   y[i] ~ normal(b[1]+b[2]*x[i]+b[3]*x[i]^2, sigma);
  // }
  y ~ normal(b[1] + b[2]*x + b[3] * x .* x, sigma);
}
'

stan_model_2 = stan_model(model_code=stan_code_2)
# save(stan_model_2, file="stan_model_2.RData")
# load("stan_model_2.RData")

fit_2  = sampling(stan_model_2,
                data=data,
                chains=3,
                iter=2000,
                warmup=1000
)

print(fit_2, probs=c(0.025, 0.975))

posterior_2 = As.mcmc.list(fit_2)

plot(posterior_2[, c("b[1]", "b[2]", "b[3]", "sigma")])

plot(posterior_2[, c("b[3]")])


# There does not seem to be much evidence for a quadratic effect, 
# there is a lot of probability mass distributed around zero 
# (more on that tomorrow).
# Reminder: Do not interpret lower order effects ("main effects") 
# independently from higher order effects!!

#------------------------------------------------------------------------------
# exercise 3 - lin.reg., different number of samples
#------------------------------------------------------------------------------

# we fit the same model to the same data, 
# using different number of posterior samples per chain (100, 1000, 10000)
# (1000 already saved in fit object)
# with larger sample size, we get a better approximation of the true posterior distribution

fit_few_samples  = sampling(stan_model,
                            data=data,
                            chains=3,
                            iter=1100,
                            warmup=1000
)

fit_many_samples  = sampling(stan_model,
                             data=data,
                             chains=3,
                             iter=11000,
                             warmup=1000
)

print(fit_few_samples)
print(fit)
print(fit_many_samples)

plot(As.mcmc.list(fit_few_samples)[,1:3])
plot(As.mcmc.list(fit_many_samples)[,1:3])

plot(fit_few_samples)
plot(fit_many_samples)

# posterior=as.matrix(fit)
# posterior_few_samples=as.matrix(fit_few_samples)
# posterior_many_samples=as.matrix(fit_many_samples)
# 
# density_1=density(posterior[, "b[2]"])
# density_few_samples=density(posterior_few_samples[, "b[2]"])
# density_many_samples=density(posterior_many_samples[, "b[2]"])
# 
# par(mfrow=c(1,1))
# 
# plot(density_1)
# lines(density_few_samples, col="red")
# lines(density_many_samples, col="blue")

#------------------------------------------------------------------------------
# exercise 4 - lin.reg., different number of observations
#------------------------------------------------------------------------------

# we generate new datasets using identical parameters, 
# but using different numbers of observations
# with larger number of observations, we get lower posterior standard deviations,
# i.e. lower uncertainty in parameter estimation.
# more information in the data (in the likelihood) means a more narrow posterior distribution. 

set.seed(123) # initiate random number generator for reproducability

n.few = 20
x.few = runif(n=n.few, min=0, max=1)
y.few = rnorm(n=n.few, mean=a+b*x.few, sd=sigma)

n.many = 500
x.many = runif(n=n.many, min=0, max=1)
y.many = rnorm(n=n.many, mean=a+b*x.many, sd=sigma)

par(mfrow=c(1,3))
plot(x.few,y.few)
plot(x,y)
plot(x.many,y.many)

data.few = list(n=n.few, 
                x=x.few,
                y=y.few)

data.many = list(n=n.many, 
                 x=x.many,
                 y=y.many)

fit_few_obs  = sampling(stan_model,
                        data=data.few,
                        chains=3,
                        iter=2000,
                        warmup=1000
)

fit_many_obs  = sampling(stan_model,
                         data=data.many,
                         chains=3,
                         iter=2000,
                         warmup=1000
)

print(fit_few_obs)
print(fit)
print(fit_many_obs)

# plot(As.mcmc.list(fit_few_obs))
# plot(As.mcmc.list(fit_many_obs))

plot(fit_few_obs)
plot(fit_many_obs)

# We extract the posterior density for the slope `b[2]` 
# and see how the posterior changes with size of the data.
# (More on handling the posterior distribution tomorrow)

posterior=as.matrix(fit)
posterior_few_obs=as.matrix(fit_few_obs)
posterior_many_obs=as.matrix(fit_many_obs)

density_1=density(posterior[, "b[2]"])
density_few_obs=density(posterior_few_obs[, "b[2]"])
density_many_obs=density(posterior_many_obs[, "b[2]"])

par(mfrow=c(1,1))

plot(density_1, xlim=c(0.5,2.5), ylim=c(0,4), main="slope b[2]")
lines(density_few_obs, col="red")
lines(density_many_obs, col="blue")
legend("topleft", bty="n", legend=c("few obs","med obs","many obs"), 
       lty=c(1,1,1), col=c("red","black","blue"))

