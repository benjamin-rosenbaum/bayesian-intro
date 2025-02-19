# Practical 4: Linear models
# by Benjamin Rosenbaum

# We learn about linear models with continuous or categorical predictors, namely 
# linear regression, ANOVA, ANCOVA
# Research questions are answered via model selection (LOO), but also with 
# comparison of posterior predictions (counterfactuals, "what-if" scenarios). 
# With categorical predictors, the emmeans package is helpful here.

rm(list=ls())
library("brms")
library("bayesplot") 
library("performance")
library("ggplot2")
library("emmeans")
library("ecostats")
library("Data4Ecologists")
try(dev.off())

setwd("~/Nextcloud/Teaching brms/Practical_04")

# 2 continuous predictors, additive --------------------------------------------

# same data as last session, but use additional predictor rainfall 
# scale predictors (to mean=0, sd=1). makes things easier espec. for interactions

data(globalPlants)
str(globalPlants)

globalPlants$z.lat = scale(globalPlants$lat)
globalPlants$z.rain = scale(globalPlants$rain)

ggplot(globalPlants, aes(z.lat, log(height))) + geom_point(alpha=0.5)
ggplot(globalPlants, aes(z.rain, log(height))) + geom_point(alpha=0.5)

# Deterministic part: mu = b0 + b1*lat + b2*rain
# Stochastic part:    logheight ~ Normal(mu,sigma)

# We use vaguely informative priors, we expect neg. relation with lat, 
# pos. with rain

fit.lm.add = brm(log(height) ~ z.lat + z.rain,
                 prior = 
                   prior(normal(-1,1), class=b, coef=z.lat) +
                   prior(normal(+1,1), class=b, coef=z.rain),
                 data = globalPlants)

# check convergence

summary(fit.lm.add, prior=TRUE)
plot(fit.lm.add)

# check predictions 
# conditional_effect will plot predictions against each predictor, 
# while the other one is held constant at its mean (here =0 because we scaled it)

plot(conditional_effects(fit.lm.add), 
     points=TRUE, 
     point_args=c(alpha=0.5),
     ask=FALSE)

# can also plot them by variable

plot(conditional_effects(fit.lm.add, effects="z.lat"), 
     points=TRUE, 
     point_args=c(alpha=0.5),
     ask=FALSE)
plot(conditional_effects(fit.lm.add, effects="z.rain"), 
     points=TRUE, 
     point_args=c(alpha=0.5),
     ask=FALSE)

# Although the model does not contain an interaction, "z.lat:z.rain" will plot 
# prediction for 3 levels of rain (mean-1sd, mean, mean+1sd) . 

plot(conditional_effects(fit.lm.add, effects="z.lat:z.rain", prob=0), 
     points=TRUE, 
     point_args=c(alpha=0.5),
     ask=FALSE)

# Lines are parallel in an additive model. 
# Can also plot vs the other predictor for 3 levels of lat

plot(conditional_effects(fit.lm.add, effects="z.rain:z.lat", prob=0), 
     points=TRUE, 
     point_args=c(alpha=0.5),
     ask=FALSE)

# model evaluation through posterior predictive checks

pp_check(fit.lm.add, ndraws=100)
pp_check(fit.lm.add, type="scatter_avg")
check_model(fit.lm.add, check=c("linearity","homogeneity","qq","normality"))

# We want to compare predictions (avg plant height) from 2 scenarios:
# mean lat&rain vs a tropical scenario (lat close to equator, high rainfall)

plot(globalPlants$z.lat, globalPlants$z.rain)
points(0,0,col="red", pch=16, cex=1.5)
points(-1.5,2,col="red", pch=16, cex=1.5)
lines(c(-1.5,0),c(2,0), col="red", lty=2)

# fitted just produces quantiles of predictions

fitted(fit.lm.add, newdata=data.frame(z.lat=-1.5, z.rain=2))
fitted(fit.lm.add, newdata=data.frame(z.lat= 0,   z.rain=0))

# need to extract full posterior predictive distributions and compute 
# distribution of difference

mu_mean     = posterior_epred(fit.lm.add, newdata=data.frame(z.lat= 0,   z.rain=0))
mu_tropical = posterior_epred(fit.lm.add, newdata=data.frame(z.lat=-1.5, z.rain=2))
hist(mu_tropical-mu_mean)
mean(mu_tropical-mu_mean)
quantile(mu_tropical-mu_mean, probs=c(0.05, 0.95))

# or use hypothesis function!

mus = fitted(fit.lm.add, 
             newdata=data.frame(z.lat =c(-1.5,0),
                                z.rain=c( 2  ,0)),
             summary=F)
mus = as.data.frame(mus)
names(mus) = c("tropical","mean")
hypothesis(mus, "tropical>mean")
hypothesis(mus, "tropical>mean") |> plot()


# 2 continuous predictors, interaction -----------------------------------------

# Does effect of rain change with latitude? --> interaction model
# Same priors for 2 main effects (slopes when other predictor is =0, here =mean )
# Vague prior on interaction, zero mean
# Mean centering makes main effects meaningful and prior choice simpler

fit.lm.int = brm(log(height) ~ z.lat * z.rain,
              prior = prior(normal(-1,1), class=b, coef=z.lat) +
                      prior(normal(+1,1), class=b, coef=z.rain) +
                      prior(normal( 0,1), class=b, coef=z.lat:z.rain),
              data = globalPlants)

# check convergence

summary(fit.lm.int, prior=TRUE)
plot(fit.lm.int)

# Now conditional_effects shows interaction for 3 levels of 2nd predictor

plot(conditional_effects(fit.lm.int, effects="z.lat:z.rain", prob=0), 
     points=TRUE, 
     point_args=c(alpha=0.5))
plot(conditional_effects(fit.lm.int, effects="z.rain:z.lat", prob=0), 
     points=TRUE, 
     point_args=c(alpha=0.5))
# --> some interaction, effect of rain (slope) becomes stronger (more important) 
#     with latitude

# posterior predictive checks

pp_check(fit.lm.int, ndraws=100)
pp_check(fit.lm.int, type="scatter_avg")
check_model(fit.lm.int, check=c("linearity","homogeneity","qq","normality"))

# inference: is interaction meaningful?

hypothesis(fit.lm.int, "z.lat:z.rain>0") 
# We observe a positive interaction of 91% certainty. 
# This is some (weak) evidence. Additionally, we perform model selection to 
# test if the interaction model produces better predictions than the additive model.

LOO(fit.lm.int, fit.lm.add)
# no, we can't say that interaction model produces better preds than additive model
# data does not sufficiently support hypothesis about interaction

# 1 categorical predictor ------------------------------------------------------

# New dataset, bird species richness in different landscapes (Data4Ecologists package)
# (1-way ANOVA)

# Question: Does species richness change with landscape type?

data(birds)
ggplot(birds, aes(x=landscape, y=S)) +
  geom_jitter(width=0.1, alpha=0.5) 

# Can either dummy-code (default) or effect-code the model (0+... removes "intercept")

default_prior(S ~ landscape, data = birds)
default_prior(S ~ 0+landscape, data = birds)
# for dummy-coding, prior for intercept (reference level mean), but not on effects!
# for effect-coding, no priors on group-level means
# --> assign meaningful priors

fit.anova1.dummy = brm(S ~ landscape,
                       prior = prior(normal(0,10), class=b),
                       data = birds)
fit.anova1.effect = brm(S ~ 0+landscape,
                        prior = prior(normal(20,10), class=b),
                        data = birds)

summary(fit.anova1.dummy, prior=TRUE)
summary(fit.anova1.effect, prior=TRUE)

# both models produce the same predictions

plot(conditional_effects(fit.anova1.dummy))
plot(conditional_effects(fit.anova1.effect))

# effect-coding gives group-level means (they are the parameters here)

summary(fit.anova1.effect)

# for dummy-coding, we can use emmeans to get these predictions

emmeans(fit.anova1.dummy, ~landscape)

# also their differences / contrasts

pairs(emmeans(fit.anova1.dummy, ~landscape))

# model selection: are there meaningful differences between landscape types?
# test against intercept-only model 

fit.anova1.null = brm(S ~ 1,
                      data = birds)

LOO(fit.anova1.null, fit.anova1.dummy)
# --> yes, strong support for ~landscape model

# make nice plots, extract ggplot object:

p1 = plot(conditional_effects(fit.anova1.dummy,"landscape"),
          points=T,
          point_args=list(width=0.1, alpha=0.3),
          plot = F
)
p1[[1]] + geom_line(aes(group="landscape"))

# model evaluation

pp_check(fit.anova1.dummy, ndraws=100)
pp_check(fit.anova1.dummy, type="scatter_avg")
check_model(fit.anova1.dummy, check=c("linearity","homogeneity","qq","normality"))

# 2 categorical predictors -----------------------------------------------------

# 2nd predictor "area", here categorical with 2 levels small / large
# Surely richness is higher in larger areas, but is the difference depending 
# on landscape type? --> fit an additive and an interaction model (dummy-coding)
# (2-way ANOVA)

birds$area = cut(birds$log.area., 2, labels=c("small","large"))
head(birds)
ggplot(birds, aes(landscape, S)) +
  geom_boxplot(aes(fill=area)) 

fit.anova2.add = brm(S ~ landscape+area,
                     prior = prior(normal(0,10), class=b),
                     data = birds)
fit.anova2.int = brm(S ~ landscape*area,
                     prior = prior(normal(0,10), class=b),
                     data = birds)

summary(fit.anova2.add, prior=TRUE)
summary(fit.anova2.int, prior=TRUE)

# Additive effects means parallel lines

p1 = plot(conditional_effects(fit.anova2.add, effects="area:landscape"),
          plot=FALSE)
p1[[1]] + geom_line(aes(group=landscape),
                    position=position_dodge(0.4))

# means and contrasts averaged over landscapes

emmeans(fit.anova2.add, ~area)
emmeans(fit.anova2.add, ~area) |> pairs()

# Additive effects means parallel lines

p1 = plot(conditional_effects(fit.anova2.add, effects="landscape:area"),
          plot=FALSE)
p1[[1]] + geom_line(aes(group=area),
                    position=position_dodge(0.4))

# means and contrasts averaged over area

emmeans(fit.anova2.add, ~landscape)
emmeans(fit.anova2.add, ~landscape) |> pairs()

# now interaction

# Interaction means individual means for all level-combinations

p1 = plot(conditional_effects(fit.anova2.int, effects="area:landscape"),
          plot=FALSE)
p1[[1]] + geom_line(aes(group=landscape),
                    position=position_dodge(0.4))

# careful when averaging over a prediction in interaction models

emmeans(fit.anova2.int, ~area)
emmeans(fit.anova2.int, ~area) |> pairs()

# Interaction means individual means for all level-combinations

p1 = plot(conditional_effects(fit.anova2.int, effects="landscape:area"),
          plot=FALSE)
p1[[1]] + geom_line(aes(group=area),
                    position=position_dodge(0.4))

# careful when averaging over a prediction in interaction models

emmeans(fit.anova2.int, ~landscape)
emmeans(fit.anova2.int, ~landscape) |> pairs()

# get means for all level combinations

emmeans(fit.anova2.int, ~landscape:area)
emmeans(fit.anova2.int, ~landscape:area) |> pairs()

LOO(fit.anova2.int, fit.anova2.add)
# no strong support for interaction model. 
# --> area-effect does not change between landscape types

# posterior predictive checks

pp_check(fit.anova2.add, ndraws=100)
pp_check(fit.anova2.add, type="scatter_avg")
check_model(fit.anova2.add, check=c("linearity","homogeneity","qq","normality"))

# categorical and continuous predictor -----------------------------------------

# now we use full resolution of area as continuous predictor. 
# does area-effect change between landscape types?
# --> test additive vs interaction ANCOVA
# as exercise: leave out priors first

ggplot(birds, aes(x=log.area., y=S, col=landscape)) +
  geom_point(alpha=0.5)  +
  facet_wrap(~landscape) +
  theme(legend.position="none")

