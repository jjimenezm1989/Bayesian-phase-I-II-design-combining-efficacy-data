get_proportion_optimal_dose_combination_above_efficacy_threshold = function(df, DLT_scenario, EFF_scenario, omega, agreement){
  
  optimal_dose_combinations_df = as.data.frame(do.call(rbind, map(df, "optimal_dose_combination")))
  
  optimal_dose_combinations_df = optimal_dose_combinations_df  %>%
    mutate(above_threshold_ind = ifelse(true_peff > p0, 1, 0)) %>%
    summarise(proportion_above_threshold = mean(above_threshold_ind)) %>%
    mutate(DLT_scenario = DLT_scenario,
           EFF_scenario = EFF_scenario,
           omega = omega,
           agreement = agreement)
  
  return(optimal_dose_combinations_df)
  
  
}