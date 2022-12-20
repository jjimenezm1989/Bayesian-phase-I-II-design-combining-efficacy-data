#' Function to calculate bias and MSE from the "df" dataset
#' @param df Dataset containing outputs (data.frame)
#' @param delta_u Design parameter (see manuscript) (numeric)
#' @param type Variable (character) to denote whether we calculate power of type-I error.
#' @param DLT_scenario Toxicity scenario (numeric)
#' @param EFF_scenario Efficacy scenario (numeric)

get_power = function(df, type, DLT_scenario, EFF_scenario, omega, agreement){
  
  delta_u = seq(0, by=0.1, 1)
  
  power = numeric()
  
  for(i in 1:length(delta_u)){
    power[i] = mean(ifelse(unlist(map(df, "max_post_prob_efficacy_higher_theta_E")) >= delta_u[i], 1, 0))
  }

  
  output_df = data.frame(DLT_scenario = DLT_scenario,
                         EFF_scenario = EFF_scenario,
                         type = type,
                         value = power,
                         delta_u = delta_u,
                         omega = omega,
                         agreement = agreement)
  
  
  output_df
  
}
