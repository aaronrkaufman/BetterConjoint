#' One of two main functions in the package.
#'
#' Given the data from a conjoint study as a vector of measurement error assumptions, conducts SIMEX to de-bias AMCE estimates.
#'
#' @param dat A data frame as per the cjoint package
#' @param formula An estiamting formula as per the cjoint package
#' @param respondent.id The name of the variable that indicates the respondent ID, as per the cjoint package
#' @param outcome.var The name of the variable that indicates respondents' choice of profile
#' @param err A vector of equal length to the number of possible conjoint profiles, representing assumptions about how much measurement error exists for each profile.
#' @param true.betas Optionally, a vector of true AMCEs
#' @param feats The vector of feature attribute lengths. Only applicable if true.betas is not null.
#' @return A SIMEX plot with the de-biased AMCE estimates
#' @export
#' @examples
#'
#'
#' @import reshape
#' @import dplyr
#' @import tidyverse
#' @import cjoint
#' @import rdist



conjoint_simex = function(dat, formula, respondent.id, outcome.var, err, true.betas=NULL, feats=NULL){
  print("Simulating error...")
  amce1<-amce(formula, data= dat, cluster=T,  respondent.id = respondent.id)
  amce2<-nperms(dat, formula, respondent.id, "choice", 1, err)
  amce3<-nperms(dat, formula, respondent.id, "choice", 2, err)
  amce4<-nperms(dat, formula, respondent.id, "choice", 3, err)
  amce5<-nperms(dat, formula, respondent.id, "choice", 4, err)

  ids = (1:length(unlist(amce1$estimates)))
  ids = ids[ids%%2 == 1]

  print("SIMEXing...")
  estimates = lapply(ids, get_one_estimate, amce1, amce2, amce3, amce4, amce5)

  if(!is.null(true.betas)){
    full.betas = unlist(lapply(1:length(true.betas), FUN=function(x) get_full_betas(feats[x], true.betas[x])))
    xmin = min(c(min(unlist(estimates), full.betas)))
    xmax = max(c(max(unlist(sapply(estimates, FUN=function(x) x[,2]))), full.betas))
  } else {
    xmin = min(unlist(estimes))
    xmax = max(unlist(sapply(estimates, FUN=function(x) x[,2])))
  }

  print("Plotting...")
  plot(NA, ylim=c(length(unlist(amce1$estimates)), 1),
       xlim=c(xmin*1.1, xmax*1.1),
       xlab="Change in E[Y]", yaxt='n', ylab="")
  axis(2, at=1:length(names(unlist(amce1$estimates))), labels=names(unlist(amce1$estimates)), las=1)
  sapply(1:length(estimates), add_points, estimates)

  if(!is.null(true.betas)){
    full.betas = unlist(lapply(1:length(true.betas), FUN=function(x) get_full_betas(feats[x], true.betas[x])))
    points(y=seq(from = length(unlist(amce1$estimates)), to = 1, by=-2)-1, x=rev(full.betas), pch=19)
  }

}
