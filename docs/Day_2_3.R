rm(list=ls())

set.seed(123)

library(rstan)
library(coda)

rstan_options(auto_write=TRUE)
options(mc.cores=4)

n=100

a=1
b=-0.5
c=-0.4
sigma=0.2

x = runif(n=n, min=-1, max=1)
y = a+b*x+c*x^2 + rnorm(n=100, mean=0, sd=sigma)

df = data.frame(x=x, y=y)

plot(df)

data = list(n=n,
            x=df$x,
            y=df$y)

stan_code = '
data {
  int n;
  vector[n] x;
  vector[n] y;
}
parameters{
  real a; // intercept
  real b; // slope
  real<lower=0> sigma; // residual sdev
}
model {
  vector[n] mu; // auxiliary variable: predictions
  // priors (weakly informative)
  a ~ normal(0,10); // mean=0, sdev=10
  b ~ normal(0,10);
  sigma ~ normal(0,5);
  
  // likelihood
  for(i in 1:n){
    mu[i] = a+b*x[i];
    y[i] ~ normal( mu[i], sigma );
  }
  // or short, in vector notation
  // mu = a+b*x;
  // y ~ normal(mu, sigma);
}
'

stan_model = stan_model( model_code=stan_code)

fit = sampling(stan_model, 
               data=data,
               iter=2000,
               warmup=1000,
               chains=3)

print(fit)

plot(fit, pars=c("a", "b"))

# stan_trace(fit)

posterior = As.mcmc.list(fit)
plot(posterior[ , 1:3])

# pairs(fit, pars=c("a","b","sigma"))

posterior = as.matrix(fit)
str(posterior)

head(posterior)

hist(posterior[ , 1])
hist(posterior[ , "b"])

pairs(fit, pars=c("a","b","sigma"))

library(BayesianTools)
correlationPlot(posterior[, 1:3], thin=1)

hist(posterior[, "b"], xlim=c(-0.7, 0))
abline(v=0)

sum(posterior[, "b"]<0) / nrow(posterior)

plot(df)

for(i in 1:50){
  abline( posterior[i,"a"], posterior[i,"b"], col=adjustcolor("red", alpha.f=0.3) )
}


x.pred=0.5
y.cred=rep(NA,50)
for(i in 1:50){
  y.cred[i] =  posterior[i,"a"] +  posterior[i,"b"]*x.pred
}

hist(y.cred)

x.pred=seq(from=-1, to=1, length=100)

x.pred
y.cred=matrix(NA, nrow=nrow(posterior), ncol=length(x.pred))
str(y.cred)
for(i in 1:nrow(posterior)){
  y.cred[i, ] =  posterior[i,"a"] +  posterior[i,"b"]*x.pred
}

y.cred[1, ]
plot(df)
lines(x.pred, y.cred[1, ])

y.cred.mean = apply(y.cred, 2, function(x) mean(x))
y.cred.q05 = apply(y.cred, 2, function(x) quantile(x, probs=0.05))
y.cred.q95 = apply(y.cred, 2, function(x) quantile(x, probs=0.95))

plot(df)
lines(x.pred, y.cred.mean, col="red")
lines(x.pred, y.cred.q05, col="red", lty=2)
lines(x.pred, y.cred.q95, col="red", lty=2)

y.pred=matrix(NA, nrow=nrow(posterior), ncol=length(x.pred))
for(i in 1:nrow(posterior)){
  y.pred[i, ] = rnorm(n=length(x.pred),
                      mean=y.cred[i, ],
                      sd=rep(posterior[i, "sigma"], length(x.pred)))
}
y.pred.mean = apply(y.pred, 2, function(x) mean(x))
y.pred.q05 = apply(y.pred, 2, function(x) quantile(x, probs=0.05))
y.pred.q95 = apply(y.pred, 2, function(x) quantile(x, probs=0.95))

lines(x.pred, y.pred.mean, col="blue")
lines(x.pred, y.pred.q05, col="blue", lty=2)
lines(x.pred, y.pred.q95, col="blue", lty=2)


x.pred = df$x
y.cred=matrix(NA, nrow=nrow(posterior), ncol=length(x.pred))
str(y.cred)
for(i in 1:nrow(posterior)){
  y.cred[i, ] =  posterior[i,"a"] +  posterior[i,"b"]*x.pred
}
y.cred.mean = apply(y.cred, 2, function(x) mean(x))

plot(y.cred.mean, df$y, xlab="predicted", ylab="observed")
abline(0,1)

plot(y.cred.mean, df$y-y.cred.mean, xlab="predicted", ylab="residual")
abline(0,0)



