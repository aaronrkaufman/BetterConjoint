#' One of two main functions in the package. 
#' 
#' Given a feature table and some assumptions about the AMCE sizes as well as a vector of measurement error assumptions, produces estimates of divergence between estimated and true AMCEs.
#'
#' @param feats An integer vector indicating the feature table. For example, c(2,3,4) indicates a 3-attribute conjoint where the first attribute has 2 levels, the second attribute has 3 levels, and the third attribute has 4 levels
#' @param min The minimum AMCE (AMCEs are assumed to derive from the uniform distribution)
#' @param max The maximum AMCE
#' @param p A vector of equal length to the number of possible conjoint profiles, representing assumptions about how much measurement error exists for each profile.
#' @return A list of length 4, each a statistic about the error in AMCE. The first element is correlation, the second is MSE, the third is MAD, the fourth is the proporiton of estimated AMCEs more extreme than the true AMCEs, and the fifth is the proportion of sign switches in the AMCE.
#' @export
#' @examples
#' feats = c(2,3,4)
#' simulate_conjoint(feats = feats, min = -0.5, max = 0.5, p=rep(0.2, prod(feats)))
#' @import reshape
#' @import dplyr
#' @import tidyverse



simulate_conjoint = function(feats, min, max, p){
  reps = replicate(100,
                   generate_data(feats=feats, min=min, max=max, p=p),
                   simplify=FALSE)
  reps2 = do.call(rbind, reps)
  out = colMeans(reps2)
  return(out)
}



