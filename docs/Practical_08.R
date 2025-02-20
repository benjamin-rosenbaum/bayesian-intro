# Practical 8: Exercise with an original dataset
# by Benjamin Rosenbaum

rm(list=ls())
library("lme4")
library("brms")
library("GGally")
library("emmeans")

setwd("~/Nextcloud/Teaching brms/Practical_08/Gould_2025_EvoEco")

# paper https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-024-02101-x#Sec2
# data source: https://osf.io/hdv8m 
# metadata: https://osf.io/hdv8m 

# In a recent study, 2 datasets were given to 174 research teams with the task to
# answer 2 research question, each. The authors collected results in form of 
# standarized effect sizes and CIs, and performed additional meta-analyses.

# We here want to "replicate" the study with the 1st dataset "evolutionary ecology".
# "To what extent is the growth of nestling blue tits (Cyanistes caeruleus) 
# influenced by competition with siblings?" Report back to me if you found a positive,
# negative, or no effect of sibling competition on individual growth. 

# Each row in the data is 1 chick individual. Note that some chicks were transferred
# to other nests after hatching, to manipulate brood size.

# - What is an appropriate response variable?
# - Which predictors measure "competition with siblings"?
# - What other predictors have to be controlled for?
# - Are there important interactions?
# - Observations are not independent: Random effects and their structure?


# Some variables I looked into (grouping by me):

# day_14_tarsus_length
# day_14_weight
# 
# hatch_Area
# hatch_Box
# hatch_mom_Ring
# hatch_nest_dad_Ring
# hatch_nest_CS
# d0_hatch_nest_brood_size
# d14_hatch_nest_brood_size 
# 
# rear_area
# rear_Box
# rear_mom_Ring
# rear_dad_Ring
# rear_nest_CS
# rear_d0_rear_nest_brood_size
# d14_rear_nest_brood_size
# 
# rear_nest_trt                 5: increase 6: decrease 7: no net brood manipulation
# home_or_away                  1: home (not transferred) 2: away (was transferred)
# net_rearing_manipulation
# day14_measurer
# chick_sex_molec
# hatch_year

# use chains=3, cores=3, to run chains in parallel for speedup
# hint: because some LMMs can take ~10 min even with parallelization,
# run fixed effects models first and then add a random effects structure

data = read.csv("blue_tit_data_updated_2020-04-18.csv", 
                na.strings=".", 
                stringsAsFactors=TRUE)

str(data)
# transform some integers to factors
data$hatch_year = as.factor(data$hatch_year)
data$rear_nest_trt = as.factor(data$rear_nest_trt)
data$day14_measurer = as.factor(data$day14_measurer)
data$home_or_away = as.factor(data$home_or_away)
data$chick_sex_molec = as.factor(data$chick_sex_molec)

# I used some variables only to test fixed effects first. I had to make sure that
# each observation is complete for model comparison (loo)
data.f = data[, c("day_14_tarsus_length",
                  "day_14_weight",
                  "hatch_nest_CS",
                  "d0_hatch_nest_brood_size",
                  "d14_hatch_nest_brood_size",
                  "rear_nest_CS",
                  "rear_d0_rear_nest_brood_size",
                  "d14_rear_nest_brood_size",
                  "rear_nest_trt",
                  "home_or_away",
                  "net_rearing_manipulation",
                  "chick_sex_molec"
)]
# remove NA observations
data.f = data.f[complete.cases(data.f), ]

ggpairs(data.f[, c("hatch_nest_CS",
                   "d0_hatch_nest_brood_size",
                   "d14_hatch_nest_brood_size",
                   "rear_nest_CS",
                   "rear_d0_rear_nest_brood_size",
                   "d14_rear_nest_brood_size",
                   "net_rearing_manipulation",
                   "day_14_tarsus_length",
                   "day_14_weight")] )


