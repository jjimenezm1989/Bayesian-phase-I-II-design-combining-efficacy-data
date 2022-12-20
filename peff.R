#'Function that calculates marginal probability of efficacy
#'
#' @param beta0 Parameter beta0 from marginal efficacy model
#' @param beta1 Parameter beta1 from marginal efficacy model
#' @param beta2 Parameter beta2 from marginal efficacy model
#' @param beta3 Parameter beta3 from marginal efficacy model


peff = function(beta0,beta1,beta2,beta3,x,y){
  p = 1/(1 + exp(-beta0 - exp(beta1)*x - exp(beta2)*y - beta3*x*y))
  return(p)
}