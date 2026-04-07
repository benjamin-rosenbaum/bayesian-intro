# Practical 8: Special topics
# by Benjamin Rosenbaum

# Zero inflated models

rm(list=ls())
library("brms")
library("ggplot2")
library("DHARMa.helpers")
library("emmeans")
library("glmmTMB") # Salamander dataset

# setwd("~/Nextcloud/Teaching brms 2025 2/00_course_folder")

# simple example ---------------------------------------------------------------

# if your data contains many zeros, this is not necessarily zero-inflation:
# look at these poisson distributions with different mean values

y = 0:20

dens = dpois(y, lambda=8)
barplot(dens, names.arg=y)

dens = dpois(y, lambda=4)
barplot(dens, names.arg=y)

dens = dpois(y, lambda=0.5)
barplot(dens, names.arg=y)

# a zero-inflated poisson is used to describe excess zeros that are not in line
# with the rest of the data. they are also called structural zeros. 
# in this example, we count species richness of an insect species that lives only
# in forests. we simulate 100 sites in forests, and 20 in grasslands. 
# if the data wouldn't include the information "location" and we want to model 
# richness in suitable sites only, just fitting a poisson model will give us a
# biased estimate for mean richness 

set.seed(100)

df = data.frame(y = rpois(100, lambda=4),
                location = "forest")
hist(df$y, breaks = seq(from=-0.5, to=15.5, by=1))

df.ZI = rbind(df, data.frame(y = rep(0,20),
                             location = "grassland"))
hist(df.ZI$y, breaks = seq(from=-0.5, to=15.5, by=1))

# we fit a poisson, a negative binomial, and a zero-inflated poisson.
# it contains a parameter "zi" in (0,1) that estimates the number of structural zeros.
# the poisson mean (and, if you have predictors, their effect sizes), are now 
# independent of the structural zeros

# P_ZIPoisson(y=0) = zi + (1-zi)*P_Poisson(y=0)
# P_ZIPoisson(y=k) = (1-zi)*P_Poisson(y=k), k>0

default_prior(y ~ 1,
              family = poisson(link="identity"),
              data = df.ZI)
default_prior(y ~ 1,
              family = negbinomial(link="identity"),
              data = df.ZI)
default_prior(y ~ 1,
              family = zero_inflated_poisson(link="identity"),
              data = df.ZI)

fit.1 = brm(y ~ 1,
            family = poisson(link="identity"),
            data = df.ZI,
            backend = "cmdstan")

fit.2 = brm(y ~ 1,
            family = negbinomial(link="identity"),
            data = df.ZI, 
            backend = "cmdstan")

fit.3 = brm(y ~ 1,
            family = zero_inflated_poisson(link="identity"),
            data = df.ZI,
            backend = "cmdstan")

dh_check_brms(fit.1)
dh_check_brms(fit.2)
dh_check_brms(fit.3)

# distributional assumptions are not met in the poisson. negative binomial does 
# not produce significant violations of assumptions, but qq-plot does not look good.
# ZIP residuals look fine!

summary(fit.1)
summary(fit.2)
summary(fit.3)

# note how the poisson, but also the negbinomial underestimate the true mean of 4.
# ZIP is unbiased. it estimates 16% structural zeros, which is pretty close to the
# actual number 20/120.

fitted(fit.3) # predictions on response scale, includes predictions for excess zeros

# note that predictions with the ZIP make predictions on the response scale, this
# means that it predicts the mean for the whole dataset including structural zeros.
# for predictions on the scale of suitable locations only, always use dpar="mu":

fitted(fit.3, dpar="mu") # poisson part only
fitted(fit.3, dpar="zi") # proportion of excess zeros 

# salamander data --------------------------------------------------------------

#  https://doi.org/10.1111/1365-2664.12585 
# (original analysis uses N-mixture model, here only simplified GLM analysis)

# salamander counts for 7 species. different sites. test effect of mined (yes/no).
# multiple locations per site. but some locations just not suitable for salamanders
# (for whatever reason). we'd want to exclude these sites, but can't know a-priori
# if a 0 is "true 0" (generally suitable, but 0 observations because small population)
# or if it is a "structural 0".

library("glmmTMB")
data(Salamanders)
str(Salamanders)
hist(Salamanders$count)

ggplot(aes(x=spp, y=count), data=Salamanders) +
  geom_jitter(width=0.1, alpha=0.3)
# we exclude 1 potential outlier
df = subset(Salamanders, count<30)

# 1. poisson -------------------------------------------------------------------

# model counts per species, effect of mining, and random intercept for site.
# poisson uses a log link, which means all additive parts in the model are 
# multiplicative on the response scale

default_prior(count ~ spp * mined + (1|site),
              family = poisson(),
              data = df)

# define additional prior for level-specific differences to reference level

fit.1 = brm(count ~ spp * mined + (1|site),
            family = poisson(),
            prior = prior(normal(0,2), class="b"),
            cores = 4,
            backend = "cmdstan",
            data = df)
summary(fit.1)
dh_check_brms(fit.1)
# strong deviation from model assumptions

plot(conditional_effects(fit.1, effects="spp:mined"),
     points=T, point_args=c(width=0.1, alpha=0.2))
# these predictions are likely biased because of the excess zeros in unsuitable locations

# 2. ZIP -----------------------------------------------------------------------

# add zi-part. note that we specify an additional model for zi, which can also
# be omitted if fitting an intercept only

default_prior(bf(count ~ spp * mined + (1|site),
                 zi ~ 1),
              family = zero_inflated_poisson(),
              data = df)

fit.2 = brm(bf(count ~ spp * mined + (1|site),
               zi ~ 1),
            family = zero_inflated_poisson(),
            prior = prior(normal(0,2), class="b"),
            cores = 4,
            backend = "cmdstan",
            data = df)
summary(fit.2)
dh_check_brms(fit.2)
# residual checks much better now

plot(conditional_effects(fit.1, effects="spp:mined"),  
     points=T, point_args=c(alpha=0.2, width=0.1))     # old
plot(conditional_effects(fit.2, effects="spp:mined"),
     points=T, point_args=c(alpha=0.2, width=0.1))     # new but wrong
plot(conditional_effects(fit.2, effects="spp:mined", dpar="mu"),
     points=T, point_args=c(alpha=0.2, width=0.1))     # new, correct
# remember to specify dpar="mu" with all predictions. only the do we see means 
# for suitable locations only.

# summary table gave us zi on logit-scale. make precition for zi to see on response
# scale, 29% excess zeros in total
fitted(fit.2, dpar="zi")
# plot(conditional_effects(fit.2, dpar="zi"))

loo(fit.1, fit.2)

# 2. ZIP with zi-model ---------------------------------------------------------

# maybe excess zeros (non-suitable locations) varies between species?

default_prior(bf(count ~ spp*mined + (1|site),
                 zi ~ spp),
              family = zero_inflated_poisson(),
              data = df)

# add a prior for dummy-coded effects for zi-model as well (on logit link)
# could also specify "zi~0+spp" to get rid of dummy-coding

fit.3 = brm(bf(count ~ spp*mined + (1|site),
               zi ~ spp),
            family = zero_inflated_poisson(),
            prior = c(prior(normal(0,2), class="b"),
                      prior(normal(0,2), class="b", dpar="zi")),
            cores = 4,
            backend = "cmdstan",
            data = df)
summary(fit.3)
dh_check_brms(fit.3)

# plot(fit.3)

plot(conditional_effects(fit.3, effects="spp:mined", dpar="mu"),
     points=T, point_args=c(alpha=0.2, width=0.1))
plot(conditional_effects(fit.3, effect="spp", dpar="zi"))

loo(fit.1, fit.2, fit.3)

# emmeans ----

# compute species-specific effects of mining on populations in suitable locations

emmeans(fit.3, ~mined|spp, type="response", dpar="mu")
emmeans(fit.3, ~mined|spp, type="response")
# emmeans does not seem to care about dpar="mu" and always predicts in suitable
# locations

emmeans(fit.3, ~mined|spp, type="response") |> pairs() # response-ratios
emmeans(fit.3, ~mined|spp) |> pairs()                  # log-differences
emmeans(fit.3, ~mined|spp) |> regrid() |> pairs()      # response-differences

emmeans(fit.3, ~mined|spp, type="response") |> pairs() |> plot()

# 3. ZIP with zi-model and cont. predictor -------------------------------------

# can suitability be predicted by "cover"?

default_prior(bf(count ~ spp*mined + (1|site),
                 zi ~ spp+cover),
              family = zero_inflated_poisson(),
              data = df)

fit.4 = brm(bf(count ~ spp*mined + (1|site),
               zi ~ spp+cover),
            family = zero_inflated_poisson(),
            prior = c(prior(normal(0,2), class="b"),
                      prior(normal(0,2), class="b", dpar="zi")),
            cores = 4,
            backend = "cmdstan",
            data = df)
summary(fit.4)
dh_check_brms(fit.4)

plot(conditional_effects(fit.4, effect="cover:spp", dpar="zi"))
plot(conditional_effects(fit.4, effect="cover:spp", dpar="zi", prob=0))
# high uncertainty and very weak effect of cover on non-suitability, but strong 
# differences between species (as in previous model)

loo(fit.3, fit.4)
# should not include cover effect on suitability. 


