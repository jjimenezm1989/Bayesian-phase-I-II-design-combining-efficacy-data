#' Function to average posterior probabilities of early stopping for toxicity, efficay and futility from the "df" dataset
#' @param df Dataset containing outputs (data.frame)
#' @param modeling Type (character) of modeling approach (Independent modeling, Latent variables, FGM Copula, Clayton Copula) (character)
#' @param delta_Z,delta_E0,delta_E1 Design parameter (see manuscript) (numeric)
#' @param DLT_scenario Toxicity scenario (numeric)
#' @param EFF_scenario Efficacy scenario (numeric)

get_posterior_probabilities_stopping_rules = function(df, delta_Z1, delta_Z2, delta_E0, DLT_scenario, EFF_scenario, omega, agreement){
  
  #Mean posterior probability early stopping for toxicity
  
  interim_prob_dlt_st2 = as_tibble(t(as.data.frame(do.call(rbind, map(df, "interim_prob_dlt_st2")))), rownames = NA)
  
  mean_post_prob_early_stopping_toxicity_st2 = mean(ifelse(unlist(lapply(interim_prob_dlt_st2,indicator_function <- function(x) sum(ifelse(x >= delta_Z2, 1, 0))))> 0, 1, 0))
  
  #Mean posterior probability early stopping for efficacy and futility
  
  interim_prob_futility = as_tibble(t(as.data.frame(do.call(rbind, map(df, "interim_prob_futility")))), rownames = NA)
  mean_post_prob_early_stopping_futility = mean(ifelse(unlist(lapply(interim_prob_futility,indicator_function <- function(x) sum(ifelse(x < delta_E0, 1, 0))))> 0, 1, 0))
  
  output_df = data.frame(DLT_scenario = DLT_scenario,
                         EFF_scenario = EFF_scenario,
                         omega = omega,
                         agreement = agreement,
                         type = c("Toxicity stopping (St. II)","Futility stopping"),
                         probability = c(mean_post_prob_early_stopping_toxicity_st2,
                                         mean_post_prob_early_stopping_futility))
  
  output_df
}









