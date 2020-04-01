### Measurement error model helpers

## Error is uniform
err_uniform = function(pairs, min, max){
  err = runif(n = nrow(pairs), min=0.5, max=1)
  return(err)
}
## Error is a function of dispersion
## That is, how similar are the two vectors
## More similar vectors are harder to tell apart, so the coefficient should be positive

err_dispersion = function(pairs, coefficient = 0.5, intercept = 0, method = "euclidean"){
  # Set up the data
  set1 = pairs %>% select(contains("p1_"))
  set2 = pairs %>% select(contains("p2_"))
  colnames(set1) = colnames(set2) = gsub("p[12]_", "p_", colnames(set1))

  # Calculate distances
  distances = sapply(1:nrow(set1), FUN=function(x) {
    tmp = rbind(set1[x,], set2[x,])
    rdist::rdist(tmp, metric=method)
  })

  # Confine it to [0.5, 1]
  err = intercept + coeffient * distances
  err[err>1] = 1
  err[err<.5] = 0.5

  return(err)
}

## Error is a function of consistency
## That is, are the features logically coherent?
# The direction argument specifies which direction features are more coherent in
# For example, is the baseline for Feature 1 liberal, while the baseline for Feature 2 is conservative?
# It should be same length as the number of features, and either 1 or -1

err_consistency = function(pairs, coefficient = 1/8, intercept = .5, direction){
  # Set up the data
  set1 = pairs %>% select(contains("p1_"))
  set2 = pairs %>% select(contains("p2_"))

  # Calculate consistency
  set1$sum = as.matrix(set1) %*% direction
  set2$sum = as.matrix(set2) %*% direction

  set1$consistency = sapply(set1$sum, FUN = function(x) {min(c(
    abs(min(set1$sum) - x),
    abs(max(set1$sum) - x)
  ))
  })
  set1$consistency = max(set1$consistency) - set1$consistency

  set2$consistency = sapply(set2$sum, FUN = function(x) {min(c(
    abs(min(set2$sum) - x),
    abs(max(set2$sum) - x)
  ))
  })
  set2$consistency = max(set2$consistency) - set2$consistency


  pair_consistency = set1$consistency + set2$consistency

  # Confine it to [0.5, 1]
  err = intercept + coeffient * pair_consistency
  err[err>1] = 1
  err[err<.5] = 0.5

  return(err)

}

## Error is a function of complexity
## That is, how much information is there?
# Here direction specifies whether features get more or less complex they leave baseline

err_complexity = function(pairs, coefficient = 1/16, intercept = .5, direction){
  # Set up the data
  set1 = pairs %>% select(contains("p1_"))
  set2 = pairs %>% select(contains("p2_"))

  # Calculate complexity
  set1$sum = as.matrix(set1) %*% direction
  set2$sum = as.matrix(set2) %*% direction

  pair_complexity = set1$sum + set2$sum
  pair_complexity = pair_complexity - min(pair_complexity)

  # Confine it to [0.5, 1]
  err = intercept + coeffient * pair_consistency
  err[err>1] = 1
  err[err<.5] = 0.5

  return(err)
}




#### SIMEX helpers

melt_pairs = function(sample_data){
 set1 = sample_data %>% select(contains("p1_"), err, respondent.id, choice)
 set2 = sample_data %>% select(contains("p2_"), err, respondent.id, choice)
 set2$choice = as.numeric(!set2$choice)
 colnames(set1) = colnames(set2) = gsub("p[12]_", "p_", colnames(set1))
 out = rbind(set1, set2)
 return(out)
}

swap01 = function(x, p){
  idx = which(rbinom(length(x), 1, prob=p)==1)
  x = as.logical(x)
  newx = x
  newx[idx] = !x[idx]
  return(as.integer(newx))
}

nperms = function(dat, formula, respondent.id, outcome.var, n, err){
  col.id = which(colnames(dat) %in% outcome.var)
  while(n > 0){
    dat[,col.id] = swap01(x=dat[,col.id], p=1-err)
    print(n)
    n = n-1
  }
  ests = amce(formula,
              cluster=T, data=dat,  respondent.id = respondent.id)
  return(ests)
}


get_one_estimate = function(id, amce1, amce2, amce3, amce4, amce5){ # note id can only be odd
  dv = c(unlist(amce1$estimates)[id],unlist(amce2$estimates)[id],
         unlist(amce3$estimates)[id],unlist(amce4$estimates)[id],
         unlist(amce5$estimates)[id])
  iv = c(1:5)
  m = lm(dv ~ iv)
  new.est = m$coefficients[1]
  out = cbind(iv, dv)
  out = rbind(c(0, new.est), out)
  out = cbind(out, id)
  return(out)
}

add_points=function(x, estimates){
  dat = estimates[[x]]
  text(y=dat[,3], x=dat[,2], labels=dat[,1])
}


#### Simulation helpers


get_betas = function(feats, min, max){
  runif(min = min, max = max, n = length(feats))
}

# Draws the true AMCEs and generates all possible profiles + true propensities using the betas
get_all_profiles = function(feats, betas){
  attrs = lapply(feats, FUN=function(x) 1:x)
  grid_pair <- expand.grid(attrs)
  grid_pair$combination = 1:nrow(grid_pair)
  grid_pair$xb = as.matrix(grid_pair[,1:length(feats)]) %*% as.matrix(betas)
  grid_pair = grid_pair %>% select(combination, xb, everything())
  return(grid_pair)
}

# for testing
#profiles = get_all_profiles(feats, min, max)

get_all_pairs = function(profiles){
  all_possible_pairs <- reshape::expand.grid.df(profiles, profiles) %>%
                                 set_names(c("combination1",
                                                   "xb1",
                                                   paste0("p1_x", 1:(ncol(profiles)-2)),
                                                   "combination2",
                                                   "xb2",
                                                   paste0("p2_x", 1:(ncol(profiles)-2)))) %>%
    # remove ties (i.e., 1024 cases of the same profile within a pair)
    filter(xb1 != xb2)
  return(all_possible_pairs)
}

#all_pairs = get_all_pairs(profiles)

# Make TRUE "population" data (note: excluding ties)
# Response is based on the latent xb variable

get_population = function(all_pairs){

  idx1 = c(2,which(grepl("p1", colnames(all_pairs))))
  idx2 = c(which(grepl("xb2", colnames(all_pairs))), which(grepl("p2", colnames(all_pairs))))

  population <- bind_rows(
    all_pairs %>%
      select(idx1) %>%
      set_names("xb", paste0("x", 1:(length(idx1)-1))) %>%
      mutate(profile = 1,
             id = row_number()),
    all_pairs %>%
      select(idx2) %>%
      set_names("xb", paste0("x", 1:(length(idx1)-1))) %>%
      mutate(profile = 2,
             id = row_number()),
  ) %>%
    group_by(id) %>%
    mutate(max_xb = max(xb)) %>%
    ungroup() %>%
    mutate(response = ifelse(xb == max_xb, 1, 0)) %>%
    select(-xb, -max_xb)

  population[grepl("x", colnames(population))] <- lapply(population[grepl("x", colnames(population))], factor)
  return(population)
}

#population = get_population(all_pairs)


# Make TRUE AMCEs ---------------------------------------------------------

get_amces = function(population){
  AMCEs <- population %>%
    select(id, profile, response, contains("x")) %>%
    gather("attribute", "level", 4:ncol(.)) %>%
    mutate(attribute = str_replace_all(attribute, "x", "") %>% as.numeric()) %>%
    group_by(attribute, level) %>%
    summarise(mean = mean(response, na.rm = TRUE)) %>%
    group_by(attribute) %>%
    mutate(mean2 = mean - mean[1]) %>%
    filter(mean2!=0) %>%
    select(-mean)
  return(AMCEs)
}


generate_data = function(feats, min, max, p){
  ## Setup
  profiles = get_all_profiles(feats, min, max)
  all_pairs = get_all_pairs(profiles)
  population = get_population(all_pairs)
  AMCEs = get_amces(population)

  ## Incorporate error
  ids = as.logical(rbinom(length(x),1,p))
  population$response[ids] = as.numeric(!population$response[ids])

  ## Draw from the data and run a regression
  results <- NULL

  for (i in 1:100){

    # Sample 1000 pairs of profiles (with replacement)
    id <- data.frame(id = 1:max(population$id))
    sample <- sample_n(id, 1000, replace = TRUE) %>%
      select(id) %>%
      pull()

    sample_df <- population %>% filter(id %in% sample) %>% select(response, contains("x"))

    out <- lm_robust(response ~ 1 + ., data = sample_df) %>%
      tidy() %>%
      mutate(id = i)

    results <- bind_rows(results, out)

  }


  ## Prep to plot
  out1 <- results %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate) %>%
    group_by(term) %>%
    summarize(amce.hat = mean(estimate))

  out2 <- AMCEs %>%
    mutate(term = paste0("x", attribute, level),
           intercept = mean2) %>%
    ungroup() %>%
    select(term, intercept)

  out3 <- merge(out1, out2)

  stats = c(correlation = cor(out3$amce.hat, out3$intercept),
               mse = sqrt(mean((out3$amce.hat - out3$intercept)^2)),
               mad = mean(abs(out3$amce.hat - out3$intercept)),
               mean_extreme = mean(abs(out3$amce.hat) > abs(out3$intercept)),
               sign_switches = mean(sign(out3$amce.hat)==sign(out3$intercept))) # add sign switches
 return(stats)
}

get_full_betas = function(x, beta){ # x is an integer indicating the number of levels of a feature
  out = c()
  out[1] = beta
  if(x > 2){
    n=2
    while(n < x){
      out[n] = beta*n
      n = n+1
    }
  }
  return(out)
}
