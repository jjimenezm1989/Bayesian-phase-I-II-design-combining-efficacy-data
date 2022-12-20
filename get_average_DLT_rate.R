#' Function to calculate the average DLT rate and proportion of trials with DLT rate above theta_Z + 0.1
#' @param df Dataset containing outputs (data.frame)
#' @param modeling Type (character) of modeling approach (Independent modeling, Latent variables, FGM Copula, Clayton Copula) (character)
#' @param DLT_scenario Toxicity scenario (numeric)
#' @param EFF_scenario Efficacy scenario (numeric)

get_average_DLT_rate = function(df, DLT_scenario, EFF_scenario, omega, agreement){
  
  Z_by_trial_st1 = map(df,"Z.stg1")
  Z_by_trial_st2 = map(df,"Z.stg2")
  Z_by_trial_overall = map(df,"Z.overall")
  
  #Average DLT rate
  average_DLT_rate_st1 = mean(unlist(Z_by_trial_st1))
  average_DLT_rate_st2 = mean(unlist(Z_by_trial_st2))
  average_DLT_rate_overall = mean(unlist(Z_by_trial_overall))
  
  #Proportion of trials with DLT rate above theta_Z + 0.1
  DLT_rate_by_trial = function(x){
    ifelse(sum(as.numeric(x))/length(as.numeric(x)) > (theta + 0.1), 1, 0)
  }
  
  proportion_trial_with_DLT_rate_above_threshold_st1 = mean(unlist(lapply(Z_by_trial_st1,DLT_rate_by_trial)))
  proportion_trial_with_DLT_rate_above_threshold_st2 = mean(unlist(lapply(Z_by_trial_st2,DLT_rate_by_trial)))
  proportion_trial_with_DLT_rate_above_threshold_overall = mean(unlist(lapply(average_DLT_rate_overall,DLT_rate_by_trial)))
  
  #Output
  output_df = data.frame(DLT_scenario = DLT_scenario,
                         EFF_scenario = EFF_scenario,
                         DLT_rate_by_trial_st1 = average_DLT_rate_st1,
                         DLT_rate_by_trial_st2 = average_DLT_rate_st2,
                         average_DLT_rate_overall = average_DLT_rate_overall,
                         proportion_trial_with_DLT_rate_above_threshold_st1 = proportion_trial_with_DLT_rate_above_threshold_st1,
                         proportion_trial_with_DLT_rate_above_threshold_st2 = proportion_trial_with_DLT_rate_above_threshold_st2,
                         proportion_trial_with_DLT_rate_above_threshold_overall = proportion_trial_with_DLT_rate_above_threshold_overall,
                         omega = omega,
                         agreement = agreement)
  
  output_df
  
}