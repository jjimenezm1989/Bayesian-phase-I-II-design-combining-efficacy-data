get_optimal_dose_combination = function(df, DLT_scenario, EFF_scenario, omega, agreement){
  
  optimal_dose_combinations_df = as.data.frame(do.call(rbind, map(df, "optimal_dose_combination")))
  
  optimal_dose_combinations_df = optimal_dose_combinations_df  %>%
    dplyr::select(x,y) %>%
    mutate(DLT_scenario = DLT_scenario,
           EFF_scenario = EFF_scenario,
           omega = omega,
           agreement = agreement)
  
  return(optimal_dose_combinations_df)
  
}
