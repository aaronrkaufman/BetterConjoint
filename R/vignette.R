setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## For now, source and load libraries manually
source("helpers.R")
source("simulate_error.R")
source("conjoint_simex.R")
library(tidyverse)
library(cjoint)

# Eventually:
# library(BetterConjoint)


## This vignette is designed to illustrate the entire conjoint process

## 1) Declare a design

# Five binary features:
#feats = rep(2, 5)

# Five features, with 2, 2, 2, 3, and 4 levels respectively
feats = c(2,2,2,3,4)

## 2) Generate data
# min and max specify the range of possible AMCEs (drawn from a uniform)

betas = get_betas(feats = feats, min = -0.3, max = 0.3)
profiles = get_all_profiles(feats = feats, betas = betas)
pairs = get_all_pairs(profiles)

## 3) Propose a measurement error model

## Caling it error is a misnomer here: it is really % agreement, so higher is less error

# No error
#pairs$err = 1

# Error is ignorable and totally random
pairs$err = err_uniform(pairs = pairs, min=0.5, max=1)

# Error is a function of dispersion
#err = err_dispersion(n = nrow(pairs))

# Error is a function of consistency
#err = err_consistency(n = nrow(pairs))

# Error is a function of complexity
#err = err_complexity(n = nrow(pairs))

## 4) Add measurement error
pairs$preference = as.numeric(pairs$xb1 >= pairs$xb2)
pairs$choice = swap01(x = pairs$preference, p = 1-pairs$err)

mean(pairs$choice == pairs$preference) # Average intracoder reliability of ~0.75

## 5) Sample some observations

# 200 respondents, 10 pairs each
n.respondents = 200
n.pairs = 10
sample_data = mosaic::sample(x = pairs, size = n.respondents * n.pairs, replace = TRUE)
sample_data$respondent.id = rep(seq(1:n.respondents), each=n.pairs)


## 6) Run SIMEX
sample_data_melted = melt_pairs(sample_data)

conjoint_simex(dat = sample_data_melted, formula = choice ~ p_x1 + p_x2 + p_x3 + p_x4 + p_x5,
               respondent.id = "respondent.id", err = sample_data_melted$err,
               true.betas = betas, feats = feats)





## For testing:
dat = sample_data_melted
formula = choice ~ p_x1 + p_x2 + p_x3 + p_x4 + p_x5
respondent.id = "respondent.id"
err = sample_data_melted$err
outcome.var = "choice"
true.betas = betas
