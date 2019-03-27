rm(list=ls())
library("plot3D")

#------------------------------------------------------------------------------

# x=seq(-3,3, length.out=10)
# 
# p=dnorm(x, mean=0, sd=1)
# 
# png("pics/grid_example_1d.png")
# 
# plot(x,p, type="l", ylim=c(0,max(p)),
#      xlab="theta", ylab="density",
#      xaxt="n",
#      yaxt="n",
#      lwd=2
#      )
# for (i in 1:length(x)){
#   lines(c(x[i],x[i]),c(0,p[i]), lwd=2, lty=2)
# }
# 
# dev.off()

#------------------------------------------------------------------------------

# x=seq(-2,2, length.out=10)
# y=x
# xy=expand.grid(x=x,y=y)
# p = with(xy, exp(-x^2-y^2))
# 
# png("pics/grid_example_2d.png")
# 
# persp(x=x, 
#       y=y, 
#       z=matrix(p,length(x),length(x)),
#       xlab=expression(theta1),
#       ylab=expression(theta2),
#       zlab="density",
#       theta = -30, 
#       phi = 30,
#       lwd=2)
# 
# dev.off()

#------------------------------------------------------------------------------

library(rstan)
library(coda)

set.seed(123) # initiate random number generator for reproducability

n=100
a=1
b=2
sigma=1
x = runif(n=n, min=-1, max=1)
y = rnorm(n=n, mean=a+b*x, sd=sigma)
df = data.frame(x=x,
                y=y)
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
  sigma ~ uniform(0, 100);
  // likelihood
  y ~ normal(a+b*x, sigma);
}
'
data = list(n=n, 
            x=df$x, 
            y=df$y)

# set options to use parallelization to run several chains simultaneously
rstan_options(auto_write = TRUE)
options(mc.cores = 3) # number of CPU cores
stan_model = stan_model(model_code=stan_code)

# diverging
fit  = sampling(stan_model,
                data=data,
                chains=3,
                iter=1100,
                warmup=100,
                init=rep(list(list(a=100,
                                   b=100,
                                   sigma=1
                                   )),3),
                algorithm="NUTS"
)

png("test1.png", width=800, height=600)

traceplot(As.mcmc.list(fit)[, 1], lwd=2, lty=1)

dev.off()

# warmup
fit  = sampling(stan_model,
                data=data,
                chains=3,
                iter=200,
                warmup=10,
                init=rep(list(list(a=20,
                                   b=20,
                                   sigma=1
                )),3),
                algorithm="NUTS"
)

png("test2.png", width=800, height=600)

traceplot(As.mcmc.list(fit)[, 1], lwd=2, lty=1)

dev.off()


