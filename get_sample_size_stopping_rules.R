get_sample_size_stopping_rules = function(df, DLT_scenario, EFF_scenario, omega, agreement){
  mean_sample_size_stopping_safety = mean(unlist(map(df, "sample_size_stopping_safety_st2")))
  mean_sample_size_stopping_futility = mean(unlist(map(df, "sample_size_stopping_futility")))
  
  output = data.frame(futility = mean_sample_size_stopping_futility, 
                      safety = mean_sample_size_stopping_safety,
                      DLT_scenario = DLT_scenario,
                      EFF_scenario = EFF_scenario,
                      omega = omega,
                      agreement = agreement) 
  
  return(output)
}