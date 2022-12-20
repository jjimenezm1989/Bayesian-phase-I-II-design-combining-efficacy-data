#'Function that generates correlated binary outputs
#'
#' @param rho00 Parameter rho00 from marginal toxicity model
#' @param rho10 Parameter rho10 from marginal toxicity model
#' @param rho00 Parameter rho01 from marginal toxicity model
#' @param alpha3 Parameter alpha3 from marginal toxicity model
#' @param beta0 Parameter beta0 from marginal efficacy model
#' @param beta1 Parameter beta1 from marginal efficacy model
#' @param beta2 Parameter beta2 from marginal efficacy model
#' @param beta3 Parameter beta3 from marginal efficacy model
#' @param x Dose from compound X
#' @param y Dose from compound Y

generate_bivariate_outcome = function(rho00,rho01,rho10,alpha3,beta0,beta1,beta2,beta3,x,y){
  
  p_Z = pdlt(rho00,rho01,rho10,alpha3,x,y)
  p_E = peff(beta0,beta1,beta2,beta3,x,y)
  
  Z_output = rbinom(n = 1, size = 1, prob = p_Z)
  E_output = rbinom(n = 1, size = 1, prob = p_E)
  
  output = c(Z_output, E_output)
  
  return(output)
  
}