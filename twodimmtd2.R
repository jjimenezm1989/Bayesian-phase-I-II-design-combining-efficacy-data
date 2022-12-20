#'Function that calculates the two-dimensional MTD curve
#'
#' @param rho00 Parameter rho00 from marginal toxicity model
#' @param rho10 Parameter rho10 from marginal toxicity model
#' @param rho00 Parameter rho01 from marginal toxicity model
#' @param alpha3 Parameter alpha3 from marginal toxicity model
#' @param x Dose x

twodimmtd2 = function(rho00,rho01,rho10,alpha3,theta,x){
  
  #logistic(u) = exp(u)/(1+exp(u))
  #logistic(u)^-1 = logit(u) = log(u/(1-u))
  
  alpha0 = log(rho00/(1-rho00))
  alpha1 = log(rho10/(1-rho10)) - log(rho00/(1-rho00))
  alpha2 = log(rho01/(1-rho01)) - log(rho00/(1-rho00))
  
  a = (log(theta/(1-theta))-alpha0-alpha1*x)/(alpha2+alpha3*x)
  
  return(a)
}




