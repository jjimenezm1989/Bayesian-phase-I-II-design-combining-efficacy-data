
get_MTD_curve = function(df, DLT_scenario, EFF_scenario, omega, agreement){
  
  mean_post_rho00 = mean(unlist(map(df, "post_rho00")))
  mean_post_rho01 = mean(unlist(map(df, "post_rho01")))
  mean_post_rho10 = mean(unlist(map(df, "post_rho10")))
  mean_post_alpha3 = mean(unlist(map(df, "post_alpha3")))
  
  q0025_post_rho00 = quantile(unlist(map(df, "post_rho00")),0.025)
  q0025_post_rho01 = quantile(unlist(map(df, "post_rho01")),0.025)
  q0025_post_rho10 = quantile(unlist(map(df, "post_rho10")),0.025)
  q0025_post_alpha3 = quantile(unlist(map(df, "post_alpha3")),0.025)
  
  q050_post_rho00 = quantile(unlist(map(df, "post_rho00")),0.5)
  q050_post_rho01 = quantile(unlist(map(df, "post_rho01")),0.5)
  q050_post_rho10 = quantile(unlist(map(df, "post_rho10")),0.5)
  q050_post_alpha3 = quantile(unlist(map(df, "post_alpha3")),0.5)
  
  q0975_post_rho00 = quantile(unlist(map(df, "post_rho00")),0.975)
  q0975_post_rho01 = quantile(unlist(map(df, "post_rho01")),0.975)
  q0975_post_rho10 = quantile(unlist(map(df, "post_rho10")),0.975)
  q0975_post_alpha3 = quantile(unlist(map(df, "post_alpha3")),0.975)
  
  
  
  output_df = data.frame(DLT_scenario = DLT_scenario,
                         EFF_scenario = EFF_scenario,
                         x = seq(0,by = 0.01,1),
                         y = c(y_mean = twodimmtd2(mean_post_rho00, mean_post_rho01, mean_post_rho10, mean_post_alpha3,theta,seq(0,by = 0.01,1)),
                               y_0025 = twodimmtd2(q0025_post_rho00, q0025_post_rho01, q0025_post_rho10, q0025_post_alpha3,theta,seq(0,by = 0.01,1)),
                               y_050 = twodimmtd2(q050_post_rho00, q050_post_rho01, q050_post_rho10, q050_post_alpha3,theta,seq(0,by = 0.01,1)),
                               y_0975 = twodimmtd2(q0975_post_rho00, q0975_post_rho01, q0975_post_rho10, q0975_post_alpha3,theta,seq(0,by = 0.01,1))),
                         type = c(rep("MTD mean",101),rep("MTD 2.5% percentile",101),rep("MTD 50% percentile",101),rep("MTD 97.5% percentile",101)),
                         omega = omega,
                         agreement = agreement)

  return(output_df)
  
}