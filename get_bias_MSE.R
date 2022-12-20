#' Function to calculate bias and MSE from the "df" dataset
#' @param df Dataset containing outputs (data.frame)
#' @param modeling Type of modeling approach (Independent modeling, Latent variables, FGM Copula, Clayton Copula) (character)
#' @param DLT_scenario Toxicity scenario (numeric)
#' @param EFF_scenario Efficacy scenario (numeric)

get_bias_MSE = function(df, DLT_scenario, EFF_scenario, omega, agreement){
  
  bias_absolute_difference_rho00 = mean(unlist(map(df, "absolute_difference_rho00")))
  mse_absolute_difference_rho00 = mean(unlist(map(df, "absolute_difference_rho00"))^2)
  
  bias_absolute_difference_rho10 = mean(unlist(map(df, "absolute_difference_rho10")))
  mse_absolute_difference_rho10 = mean(unlist(map(df, "absolute_difference_rho10"))^2)
  
  bias_absolute_difference_rho01 = mean(unlist(map(df, "absolute_difference_rho01")))
  mse_absolute_difference_rho01 = mean(unlist(map(df, "absolute_difference_rho01"))^2)
  
  bias_absolute_difference_alpha3 = mean(unlist(map(df, "absolute_difference_alpha3")))
  mse_absolute_difference_alpha3 = mean(unlist(map(df, "absolute_difference_alpha3"))^2)
  
  bias_absolute_difference_beta0_st2 = mean(unlist(map(df, "absolute_difference_beta0_st2")))
  mse_absolute_difference_beta0_st2 = mean(unlist(map(df, "absolute_difference_beta0_st2"))^2)

  bias_absolute_difference_beta1_st2 = mean(unlist(map(df, "absolute_difference_beta1_st2")))
  mse_absolute_difference_beta1_st2 = mean(unlist(map(df, "absolute_difference_beta1_st2"))^2)
  
  bias_absolute_difference_beta2_st2 = mean(unlist(map(df, "absolute_difference_beta2_st2")))
  mse_absolute_difference_beta2_st2 = mean(unlist(map(df, "absolute_difference_beta2_st2"))^2)
  
  bias_absolute_difference_beta3_st2 = mean(unlist(map(df, "absolute_difference_beta3_st2")))
  mse_absolute_difference_beta3_st2 = mean(unlist(map(df, "absolute_difference_beta3_st2"))^2)
 
  output_df = data.frame(DLT_scenario = DLT_scenario,
                         EFF_scenario = EFF_scenario,
                         type = c(rep("Bias",8),rep("MSE",8)),
                         parameter = rep(c("rho00", "rho10", "rho01", "alpha3", "beta0", "beta1","beta2", "beta3"),2),
                         omega = omega,
                         agreement = agreement,
                         value = c(bias_absolute_difference_rho00,
                                   bias_absolute_difference_rho10,
                                   bias_absolute_difference_rho01,
                                   bias_absolute_difference_alpha3,
                                   bias_absolute_difference_beta0_st2,
                                   bias_absolute_difference_beta1_st2,
                                   bias_absolute_difference_beta2_st2,
                                   bias_absolute_difference_beta3_st2,
                                   mse_absolute_difference_rho00,
                                   mse_absolute_difference_rho10,
                                   mse_absolute_difference_rho01,
                                   mse_absolute_difference_alpha3,
                                   mse_absolute_difference_beta0_st2,
                                   mse_absolute_difference_beta1_st2,
                                   mse_absolute_difference_beta2_st2,
                                   mse_absolute_difference_beta3_st2))
  
  
  output_df
  
}

