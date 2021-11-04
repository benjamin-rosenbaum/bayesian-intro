rm(list=ls())

set.seed(123)

library(rstan)
library(coda)
library(BayesianTools)
library(brms)

rstan_options(auto_write=TRUE)
options(mc.cores=4)

a=1
b=0.5
c=0.4
d=0.5
sigma=0.1

n=100

x1=runif(n=n, min=0, max=2)
x2=runif(n=n, min=0, max=2)
mu = a + b*x1 + c*x2 + d*x1*x2
y = rnorm(n=n, mean=mu, sd=sigma)

df = data.frame(x1=x1,
                x2=x2,
                y=y)

head(df)

stan_code = '
data{
  int n;
  vector[n] x1;
  vector[n] x2;
  vector[n] y;
}
parameters{
  real a;
  real b;
  real c;
  real<lower=0> sigma;
}
model{
  a ~ normal(0, 10);
  b ~ normal(0, 10);
  c ~ normal(0, 10);
  sigma ~ normal(0, 10);
  // likelihood
  y ~ normal( a + b*x1 + c*x2, sigma );
}
'

data = list(n=nrow(df),
            x1 = df$x1,
            x2 = df$x2,
            y = df$y )

stan_model = stan_model(model_code = stan_code)
fit.1 = sampling(stan_model, data=data)

print(fit.1)

post = as.matrix(fit.1)
head(post)

y.fit = matrix(NA, nrow=nrow(post), ncol=data$n)
for(i in 1:nrow(post)){
  y.fit[i, ] = post[i, "a"] + post[i, "b"]*data$x1 + post[i, "c"]*data$x2
}
head(y.fit)

y.fit.mean=colMeans(y.fit)

plot(y.fit.mean, df$y, xlab="pred", ylab="obs")
abline(0,1)

plot(y.fit.mean, df$y-y.fit.mean, xlab="pred", ylab="res")
abline(0,0)

x1.pred = seq(from=min(df$x1), to=max(df$x1), length.out=100)
x2.pred = mean(df$x2)
y.fit = matrix(NA, nrow=nrow(post), ncol=100)
for(i in 1:nrow(post)){
  y.fit[i, ] = post[i, "a"] + post[i, "b"]*x1.pred + post[i, "c"]*x2.pred
}
y.fit.mean=colMeans(y.fit)

plot(df$x1, df$y)
lines(x1.pred, y.fit.mean, col="red")

x2.pred = mean(df$x2)-sd(df$x2)
y.fit = matrix(NA, nrow=nrow(post), ncol=100)
for(i in 1:nrow(post)){
  y.fit[i, ] = post[i, "a"] + post[i, "b"]*x1.pred + post[i, "c"]*x2.pred
}
y.fit.mean=colMeans(y.fit)
lines(x1.pred, y.fit.mean, col="green")





