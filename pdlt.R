#'Function that calculates marginal probability of DLT
#'
#' @param rho00 Parameter rho00 from marginal toxicity model
#' @param rho10 Parameter rho10 from marginal toxicity model
#' @param rho00 Parameter rho01 from marginal toxicity model
#' @param alpha3 Parameter alpha3 from marginal toxicity model


pdlt = function(rho00,rho01,rho10,alpha3,x,y){
  alpha0 = log(rho00/(1-rho00))
  alpha1 = log(rho10/(1-rho10)) - log(rho00/(1-rho00))
  alpha2 = log(rho01/(1-rho01)) - log(rho00/(1-rho00))
  p = 1/(1+exp(-alpha0-alpha1*x-alpha2*y-alpha3*x*y))
  return(p)
}