# Practical 6: Mixed effects models
# by Benjamin Rosenbaum

# We perform some classical linear modeling with a random grouping factor. 
# Predictions etc can be computed with random effects, or fixed effects only.
# Model checks, model comparisons, hypotheses, etc are the same as with linear models.

rm(list=ls())
library("brms")
library("ggplot2")
library("performance")
library("ecostats")

setwd("~/Nextcloud/Teaching brms/Practical_06")


# random int. ANOVA ------------------------------------------------------------

# from Warton, D. (2022) Eco-Stats: Data Analysis in Ecology

# number of invertebrate in a sample from different estuaries.
# estuaries are either pristine or modified.
# test for difference between mod/prist while account for grouping factor Estuary
# "Estuary" actually nested in "Mod" (Mod/Estuary)
# Here we can treat Mod as fixed and Estuary as random effect

data("estuaries")
table(estuaries$Mod, estuaries$Estuary)
ggplot(estuaries, aes(Estuary, Total, col=Mod)) + geom_jitter(width=0.1, size=2)

# Remove NAs. brms automatically removes NA observations (rows in dataframe).
# But when they appear in different predictors, make sure that all models are 
# fitted on the same dataset. Otherwise, model comparison doesn't work.

ID.complete = complete.cases(estuaries[, c("Total","Mod","Estuary","Temperature")])
estuaries = estuaries[ID.complete, ]

# partial pooling model

# deterministic part: mu[i] = b[Mod[i]] + delta[Estuary[i]]   i=1:N
#                     Mod = 1:2, Estuary=1:7
# stochastic part:    Total[i] ~ normal(mu[i], sigma)         i=1:N
# hierarchical:       delta[j] ~ normal(0, sd_Estuary)        j=1:7
# priors:             b ~..., sd_Estuary~..., sigma~...

# (here b is effect-coded, but the fitted model below is dummy-coded)

default_prior(Total ~ Mod+(1|Estuary), data=estuaries)

fit.est.1 = brm(Total ~ Mod+(1|Estuary),
                prior = prior(normal(0,10), class=b),
                data = estuaries)
# some divergent transitions, but just 3 out of 4000 draws is no big deal
summary(fit.est.1, prior=TRUE)
plot(fit.est.1)

# extract fixed and random effects (random effects not listed in summary)

fixef(fit.est.1)
ranef(fit.est.1)

# fixed effects predictions

plot(conditional_effects(fit.est.1), 
     points=T, 
     point_args=c(alpha=0.3, width=0.1))

# fixed & random effects predictions
# We get also get predictions (on Estuary-level) of what would happen if we 
# applied modification instead of keeping an estuary pristine.

plot(conditional_effects(fit.est.1,
                         re_formula = NULL,
                         conditions = make_conditions(fit.est.1, var=c("Estuary"))),
     points=T, 
     point_args=c(alpha=0.3, width=0.1, size=2))

# when computing predictions, we can choose to include random effects (default),
# or fixed effects only (re_formula=NA)
# re_formula=NULL **enforces** random effects, re_formula=NA **omits** random effects.

fitted(fit.est.1) |> head()
fitted(fit.est.1, re_formula=NA) |> head()

# The same holds for, e.g., R2 values. With random effects (default) computes the 
# conditional R2, while without random effects (re_formula=NA) computes the 
# marginal R2.

bayes_R2(fit.est.1)
bayes_R2(fit.est.1, re_formula=NULL) # (the same)
bayes_R2(fit.est.1, re_formula=NA)

# Model checks use random effects per default, but can also be done with fixed effects
# only if required.

pp_check(fit.est.1, ndraws=100)
pp_check(fit.est.1, type="scatter_avg")

pp_check(fit.est.1, ndraws=100, re_formula=NA)
pp_check(fit.est.1, type="scatter_avg", re_formula=NA)

# since this is a linear model, we can assess model assumptions with check_model
# there is an additional check for the random effects (reqq)

check_model(fit.est.1, check=c("linearity","homogeneity","qq","normality"))
check_model(fit.est.1, check=c("reqq"))

# test for difference Mod-Pristine

hypothesis(fit.est.1, "ModPristine<0", alpha=0.05)

# random int. ANCOVA ------------------------------------------------------------

# deterministic part: mu[i] = b[Mod[i]] + c*Temperature[i] + delta[Estuary[i]]  
#                     Mod = 1:2, Estuary=1:7
# stochastic part:    Total[i] ~ normal(mu[i], sigma)         i=1:N
# hierarchical:       delta[j] ~ normal(0, sd_Estuary)        j=1:7
# priors:             b ~..., c~..., sd_Estuary~..., sigma~...

# (here b is effect-coded, but the fitted model below is dummy-coded)

default_prior(Total ~ Mod+scale(Temperature)+(1|Estuary),
              data = estuaries)

# weak prior for categorical and continuous variables

fit.est.2 = brm(Total ~ Mod+scale(Temperature)+(1|Estuary),
                prior = prior(normal(0,10), class=b),
                data = estuaries)

summary(fit.est.2, prior=TRUE)
plot(fit.est.2)

# fixed effects predictions

plot(conditional_effects(fit.est.2, effect="Temperature:Mod"), 
     points=T)

# fixed & random effects predictions

plot(conditional_effects(fit.est.2, effect="Temperature:Mod",
                         re_formula = NULL,
                         conditions = make_conditions(fit.est.2, var=c("Estuary")),
                         prob=0),
     points=T, 
     point_args=c(alpha=1, width=0.1, size=2))

pp_check(fit.est.2, ndraws=100)
pp_check(fit.est.2, type="scatter_avg")

# There is very weak support for the temperature model.
# Both models come to the similar conclusions regarding Mod-Pristine

LOO(fit.est.1, fit.est.2)
hypothesis(fit.est.1, "ModPristine<0")
hypothesis(fit.est.2, "ModPristine<0")

# exercise ---------------------------------------------------------------------

# from Warton, D. (2022) Eco-Stats: Data Analysis in Ecology

# test additional effect of outer vs inner (Zone) in new dataset
# does it hold when accounting for temperature ?

data("estuaryZone")
table(estuaryZone$Estuary, estuaryZone$Mod, estuaryZone$Zone)

# Remove NAs
ID.complete = complete.cases(estuaryZone[, c("Total","Mod","Estuary","Temperature","Zone")])
estuaryZone = estuaryZone[ID.complete, ]

