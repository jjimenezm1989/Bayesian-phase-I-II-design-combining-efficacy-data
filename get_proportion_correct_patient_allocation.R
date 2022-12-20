#' Function to calculate proportion of correct patient allocation in phase II from the "df" dataset
#' @param df Dataset containing outputs (data.frame)
#' @param modeling Type of modeling approach (Independent modeling, Latent variables, FGM Copula, Clayton Copula) (character)
#' @param DLT_scenario Toxicity scenario (numeric)
#' @param EFF_scenario Efficacy scenario (numeric)

get_proportion_correct_patient_allocation = function(df, DLT_scenario, EFF_scenario, omega, agreement){
  
  proportion = mean(unlist(map(df, "proportion_allocation_efficacious_doses")))
  
  output_df = data.frame(DLT_scenario = DLT_scenario,
                         EFF_scenario = EFF_scenario,
                         value = proportion,
                         omega = omega,
                         agreement = agreement)
  
  output_df
  
}