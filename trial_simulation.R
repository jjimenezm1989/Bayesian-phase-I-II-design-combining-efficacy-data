

trial_simulation = function(true_rho00, true_rho01, true_rho10, true_alpha3, 
                            true_beta0_st1, true_beta1_st1, true_beta2_st1, true_beta3_st1,
                            true_beta0_st2, true_beta1_st2, true_beta2_st2, true_beta3_st2,
                            path_JAGS, w){

    
  jags_txt_stage1 = paste0(path_JAGS,"/stage1.bug.txt")
  jags_txt_stage2 = paste0(path_JAGS,"/stage2.bug.txt")
  
  
  # Declaration of the output statistics
  interim_prob_dlt_st1 = numeric()
  
  # Beginning of study simulation
  
  #Stage I:
  
  #First 2 patients
  
  doseX1 = c(0,0)
  doseY1 = c(0,0)
  sample_size_stopping_safety_st1 = m1
  stopping_safety_ind = 0
  
  
  outcomes_1 = generate_bivariate_outcome(rho00 = true_rho00, rho01 = true_rho01, 
                                          rho10 = true_rho10, alpha3 = true_alpha3,
                                          beta0 = true_beta0_st1, beta1 = true_beta1_st1, 
                                          beta2 = true_beta2_st1, beta3 = true_beta3_st1, 
                                          x = doseX1[1], y = doseY1[1])
  
  outcomes_2 = generate_bivariate_outcome(rho00 = true_rho00, rho01 = true_rho01, 
                                          rho10 = true_rho10, alpha3 = true_alpha3,
                                          beta0 = true_beta0_st1, beta1 = true_beta1_st1, 
                                          beta2 = true_beta2_st1, beta3 = true_beta3_st1, 
                                          x = doseX1[2], y = doseY1[2])
  
  Z = c(outcomes_1[1],outcomes_2[1])
  NEff.stg1 = c(outcomes_1[2],outcomes_2[2])
  Ncomb1 = m1 
  
  data_jags = list('Z'= Z, 'X'= doseX1, 'Y'= doseY1, 'theta'= theta, 'N'= Ncomb1, 'a00' = a00, 'b00' = b00, 
                   'a01' = a01, 'b01' = b01, 'a10' = a10, 'b10' = b10, 'a' = a, 'b' = b)
  jags <- jags.model(jags_txt_stage1,data = data_jags, n.chains=chains,n.adapt=burn,quiet = TRUE)
  s=coda.samples(jags,c('mtdx1','mtdx2','mtdy1','mtdy2',"rho00","rho10","rho01","alpha3"),mm,progress.bar = "none")
  ss=as.data.frame(s[[1]])
  
  #Calculating posterior probability of toxicity at x = 0, y = 0.
  #We create, for each MCMC sample, the correspondent P(DLT). 
  #We then calculate the mean posterior P(DLT) at x = 0, y = 0
  ss2 = ss %>% mutate(temp1 = pdlt(rho00,rho01,rho10,alpha3, x = 0, y = 0),
                      temp2 = ifelse(temp1 > theta + 0.1, 1, 0))
  
  interim_prob_dlt_st1[1] = dplyr::pull(ss2 %>% summarise(m = mean(temp2)))
  
  #We assess whether we will enroll a new cohort or not based on safety stopping rule
  stopping_safety_ind = ifelse(stopping_safety_ind == 0 & (interim_prob_dlt_st1[1] <= delta_Z1), 0, 1)
  sample_size_stopping_safety_st1 = ifelse(stopping_safety_ind == 0 & (interim_prob_dlt_st1[1] <= delta_Z1), sample_size_stopping_safety_st1 + 2, sample_size_stopping_safety_st1)
  
  
  trcmtdx1 = ss$mtdx1[ss$mtdx1 > lbx]
  xx1 = max(0,quantile(trcmtdx1,alpha[1]))
  xx1 = min(xx1,1)
  if ((xx1 - doseX1[1]) > delta1)
    xx1 = doseX1[1]+delta1
  
  trcmtdy2 = ss$mtdy2[ss$mtdy2 > lby]
  yy2 = max(0,quantile(trcmtdy2,alpha[1]))
  yy2 = min(yy2,1)
  if ((yy2 - doseY1[2]) > delta1)
    yy2 = doseY1[2]+delta1
  
  # Update data
  doseX1 = c(doseX1,0,xx1)
  doseY1 = c(doseY1,yy2,0)
  
  #Patients 3 and 4
  
  outcomes_1 = generate_bivariate_outcome(rho00 = true_rho00, rho01 = true_rho01, 
                                          rho10 = true_rho10, alpha3 = true_alpha3,
                                          beta0 = true_beta0_st1, beta1 = true_beta1_st1, 
                                          beta2 = true_beta2_st1, beta3 = true_beta3_st1, 
                                          x = doseX1[3], y = doseY1[3])
  
  outcomes_2 = generate_bivariate_outcome(rho00 = true_rho00, rho01 = true_rho01, 
                                          rho10 = true_rho10, alpha3 = true_alpha3,
                                          beta0 = true_beta0_st1, beta1 = true_beta1_st1, 
                                          beta2 = true_beta2_st1, beta3 = true_beta3_st1, 
                                          x = doseX1[4], y = doseY1[4])
  
  Z = c(Z,outcomes_1[1],outcomes_2[1])
  NEff.stg1 = c(NEff.stg1,outcomes_1[2],outcomes_2[2])
  Ncomb1 = Ncomb1+m1
  
  
  for (c1 in 3:(NN_1/2)){
    
    data_jags = list('Z'= Z, 'X'= doseX1, 'Y'= doseY1, 'theta'= theta, 'N'= Ncomb1, 'a00' = a00, 'b00' = b00, 
                     'a01' = a01, 'b01' = b01, 'a10' = a10, 'b10' = b10, 'a' = a, 'b' = b)
    jags <- jags.model(jags_txt_stage1, data = data_jags, n.chains=chains,n.adapt=burn,quiet = TRUE)
    s=coda.samples(jags,c('mtdx1','mtdx2','mtdy1','mtdy2',"rho00","rho10","rho01","alpha3"),mm,progress.bar = "none")
    ss=as.data.frame(s[[1]])
    
    #Calculating posterior probability of toxicity at x = 0, y = 0.
    #We create, for each MCMC sample, the correspondent P(DLT). 
    #We then calculate the mean posterior P(DLT) at x = 0, y = 0
    ss2 = ss %>% mutate(temp1 = pdlt(rho00,rho01,rho10,alpha3, x = 0, y = 0),
             temp2 = ifelse(temp1 > theta + 0.1, 1, 0))
    
    interim_prob_dlt_st1[c1-1] = dplyr::pull(ss2 %>% summarise(m = mean(temp2)))
    
    #We assess whether we will enrol a new cohort or not based on safety stopping rule
    stopping_safety_ind = ifelse(stopping_safety_ind == 0 & (interim_prob_dlt_st1[c1-1] <= delta_Z1), 0, 1)
    sample_size_stopping_safety_st1 = ifelse(stopping_safety_ind == 0 & (interim_prob_dlt_st1[c1-1] <= delta_Z1), sample_size_stopping_safety_st1 + 2, sample_size_stopping_safety_st1)
    
    
    trcmtdx1 = ss$mtdx1[ss$mtdx1 > lbx]
    xx1 = max(0,quantile(trcmtdx1,alpha[c1-1]))
    xx1 = min(xx1,1)
    if ((xx1 - doseX1[2*c1-3]) > delta1)
      xx1 = doseX1[2*c1-3]+delta1
    
    trcmtdx2 = ss$mtdx2[ss$mtdx2 > lbx]
    xx2 = max(0,quantile(trcmtdx2,alpha[c1-1]))
    xx2 = min(xx2,1)
    if ((xx2 - doseX1[2*c1-2]) > delta1)
      xx2 = doseX1[2*c1-2]+delta1
    
    trcmtdy1 = ss$mtdy1[ss$mtdy1 > lby]
    yy1 = max(0,quantile(trcmtdy1,alpha[c1-1]))
    yy1 = min(yy1,1)
    if ((yy1 - doseY1[2*c1-3]) > delta1)
      yy1 = doseY1[2*c1-3]+delta1
    
    trcmtdy2 = ss$mtdy2[ss$mtdy2 > lby]
    yy2 = max(0,quantile(trcmtdy2,alpha[c1-1]))
    yy2 = min(yy2,1)
    if ((yy2 - doseY1[2*c1-2]) > delta1)
      yy2 = doseY1[2*c1-2]+delta1
    
    
    if (doseX1[2*c1-3] == doseX1[2*c1-5]) {
      
      doseX1 = c(doseX1,xx1,doseX1[2*c1-2])
      doseY1 = c(doseY1,doseY1[2*c1-3],yy2)} else {
        doseX1 = c(doseX1,doseX1[2*c1-3],xx2)
        doseY1 = c(doseY1,yy1,doseY1[2*c1-2])}
    
    
    outcomes_1 = generate_bivariate_outcome(rho00 = true_rho00, rho01 = true_rho01, 
                                            rho10 = true_rho10, alpha3 = true_alpha3,
                                            beta0 = true_beta0_st1, beta1 = true_beta1_st1, 
                                            beta2 = true_beta2_st1, beta3 = true_beta3_st1, 
                                            x = doseX1[2*c1-1], y = doseY1[2*c1-1])
    
    outcomes_2 = generate_bivariate_outcome(rho00 = true_rho00, rho01 = true_rho01, 
                                            rho10 = true_rho10, alpha3 = true_alpha3,
                                            beta0 = true_beta0_st1, beta1 = true_beta1_st1, 
                                            beta2 = true_beta2_st1, beta3 = true_beta3_st1, 
                                            x = doseX1[2*c1-2], y = doseY1[2*c1-2])
    
    Z = c(Z,outcomes_1[1],outcomes_2[1])
    NEff.stg1 = c(NEff.stg1,outcomes_1[2],outcomes_2[2])
    Ncomb1 = Ncomb1+m1
    
  }
  
  #Fix MTD curve model parameters
  
  data_jags = list('Z'= Z, 'X'= doseX1, 'Y'= doseY1, 'theta'= theta, 'N'= Ncomb1, 'a00' = a00, 'b00' = b00, 
                   'a01' = a01, 'b01' = b01, 'a10' = a10, 'b10' = b10, 'a' = a, 'b' = b)
  jags <- jags.model(jags_txt_stage1, data = data_jags, n.chains=chains,n.adapt=burn,quiet = TRUE)
  s=coda.samples(jags,c('mtdx1','mtdx2','mtdy1','mtdy2',"rho00","rho10","rho01","alpha3"),mm,progress.bar = "none")
  ss=as.data.frame(s[[1]])
  
  #Calculating posterior probability of toxicity at x = 0, y = 0.
  #We create, for each MCMC sample, the correspondent P(DLT). 
  #We then calculate the mean posterior P(DLT) at x = 0, y = 0
  ss2 = ss %>% mutate(temp1 = pdlt(rho00,rho01,rho10,alpha3, x = 0, y = 0),
                      temp2 = ifelse(temp1 > theta + 0.1, 1, 0))
  
  interim_prob_dlt_st1[NN_1/2] = dplyr::pull(ss2 %>% summarise(m = mean(temp2)))
  
  #We assess whether we will enroll a new cohort or not based on safety stopping rule
  stopping_safety_ind = ifelse(stopping_safety_ind == 0 & (interim_prob_dlt_st1[NN_1/2] <= delta_Z1), 0, 1)
  sample_size_stopping_safety_st1 = ifelse(stopping_safety_ind == 0 & (interim_prob_dlt_st1[NN_1/2] <= delta_Z1), sample_size_stopping_safety_st1 + m2, sample_size_stopping_safety_st1)
  
  
  post_rho00 = median(ss$rho00)
  post_rho01 = median(ss$rho01)
  post_rho10 = median(ss$rho10)
  post_alpha3 = median(ss$alpha3)
  
   
  #STAGE II
  
  wMix = c(w,(1-w))
  
  stopping_futility_ind = 0
  interim_prob_futility = numeric()
  sample_size_stopping_futility = NN_1 + m2
  
  stopping_safety_ind_st2 = 0
  interim_prob_dlt_st2 = numeric()
  sample_size_stopping_safety_st2 = n2
  
  first_cohort_st2 = data.frame(xx1 = seq(0,1,length.out = n2)) %>%
    mutate(yy1 = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3,theta,xx1),
           yy1 = ifelse(yy1 < 0, 0, ifelse(yy1 > 1, 1, yy1)))
  
  doseX2 =  first_cohort_st2$xx1
  doseY2 =  first_cohort_st2$yy1
  Ncomb2 = n2
  
  outcomes_n2 = mapply(generate_bivariate_outcome,
                       rho00 = true_rho00, rho01 = true_rho01, 
                       rho10 = true_rho10, alpha3 = true_alpha3,
                       beta0 = true_beta0_st2, beta1 = true_beta1_st2, 
                       beta2 = true_beta2_st2, beta3 = true_beta3_st2,
                       x = doseX2, y = doseY2)
                    
  Z.stg2 = outcomes_n2[1,]
  NEff.stg2 = outcomes_n2[2,]
  n_dlts = sum(Z.stg2)
  n_doses = length(Z.stg2)
  
  data_jags = list('NEff.stg1' = NEff.stg1, 'doseX1' = doseX1, 'doseY1' = doseY1, 'Ncomb1' = Ncomb1,
                   'NEff.stg2' = NEff.stg2, 'doseX2' = doseX2, 'doseY2' = doseY2, 'Ncomb2' = Ncomb2, 
                   'wMix' = wMix, 'theta' = theta, 'n_dlts' = n_dlts, 'n_doses' = n_doses)
  
  j=jags.model(jags_txt_stage2,data= data_jags, n.chains=chains,n.adapt=burn,quiet = TRUE)
  s=coda.samples(j,c('beta02','beta12','beta22','beta32','prob.ex','safety_stopping'),mm,progress.bar="none")
  ss=as.data.frame(s[[1]])
  
  post_beta0_st2 = median(ss$beta02)
  post_beta1_st2 = median(ss$beta12)
  post_beta2_st2 = median(ss$beta22)
  post_beta3_st2 = median(ss$beta32)
  post_prob.ex = mean(ss$`prob.ex[1]`)
  
  #Safety stopping rule
  interim_prob_dlt_st2[1] = dplyr::pull(ss %>% summarise(stopping_rule_stage2 = mean(safety_stopping)))
  
  #We assess whether we will enroll a new cohort or not based on safety stopping rule
  stopping_safety_ind_st2 = ifelse(stopping_safety_ind_st2 == 0 & (interim_prob_dlt_st2[1] <= delta_Z2), 0, 1)
  sample_size_stopping_safety_st2 = ifelse(stopping_safety_ind_st2 == 0 & interim_prob_dlt_st2[1] <= delta_Z2 & (sample_size_stopping_safety_st2 < (NN_1 + NN_2)), sample_size_stopping_safety_st2 + m2, sample_size_stopping_safety_st2)
  
  
  for(c2 in 1:((NN_2-n2)/m2)){
    
    integrate_peff = function(x){
      y_mtd = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3,theta,x)
      y_mtd = ifelse(y_mtd < 0, 0, ifelse(y_mtd > 1, 1, y_mtd))
      p = peff(post_beta0_st2,post_beta1_st2,post_beta2_st2,post_beta3_st2,x,y_mtd)
      return(p)
      
    }
    
    ncst<-integrate(integrate_peff,0,1)$value
    
    tempmax = dplyr::pull(data.frame(x_grid = x_grid) %>% 
                            mutate(y = integrate_peff(x_grid)) %>%
                            summarise(tempmax = max(y)))
    
    M<-(1/ncst)*tempmax
    
    #Sample from the standardized density based on the rejection-sampling principle.
    
    rejection_sampling_df = data.frame(sample.x = runif(600,0,1), U = runif(600,0,1)) %>%
      mutate(rejection_boundary = (1/ncst)*integrate_peff(sample.x),
             density_sample_value = dunif(sample.x, 0, 1)*M*U,
             acceptance = ifelse(density_sample_value <= rejection_boundary, "Yes", "No")) %>%
      filter(acceptance == "Yes") %>% 
      dplyr::select(sample.x) %>% 
      slice(1:m2) %>% 
      mutate(y = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3,theta,sample.x),
             y = ifelse(y < 0, 0, ifelse(y > 1, 1, y)))

    xx1 = rejection_sampling_df$sample.x
    yy1 = rejection_sampling_df$y

    outcomes = mapply(generate_bivariate_outcome, 
                      rho00 = true_rho00, rho01 = true_rho01, 
                      rho10 = true_rho10, alpha3 = true_alpha3,
                      beta0 = true_beta0_st2, beta1 = true_beta1_st2, 
                      beta2 = true_beta2_st2, beta3 = true_beta3_st2,
                      x = xx1, y = yy1)
    
    Z.stg2 = c(Z.stg2,outcomes[1,])
    NEff.stg2 = c(NEff.stg2,outcomes[2,])
    doseX2 = c(doseX2,xx1)
    doseY2 = c(doseY2,yy1)
    
    Ncomb2 = Ncomb2 + m2
    
    n_dlts = sum(Z.stg2)
    n_doses = length(Z.stg2)

      
    data_jags = list('NEff.stg1' = NEff.stg1, 'doseX1' = doseX1, 'doseY1' = doseY1, 'Ncomb1' = Ncomb1,
                     'NEff.stg2' = NEff.stg2, 'doseX2' = doseX2, 'doseY2' = doseY2, 'Ncomb2' = Ncomb2, 
                     'wMix' = wMix, 'theta' = theta, 'n_dlts' = n_dlts, 'n_doses' = n_doses)
    
    j=jags.model(jags_txt_stage2,data= data_jags, n.chains=chains,n.adapt=burn,quiet = TRUE)
    s=coda.samples(j,c('beta02','beta12','beta22','beta32','prob.ex','safety_stopping'),mm,progress.bar="none")
    ss=as.data.frame(s[[1]])
    
    #Calculating posterior probability of toxicity in Stage II.
    interim_prob_dlt_st2[c2 + 1] = dplyr::pull(ss %>% summarise(stopping_rule_stage2 = mean(safety_stopping)))
    
    #We assess whether we will enroll a new cohort or not based on safety stopping rule
    stopping_safety_ind_st2 = ifelse(stopping_safety_ind_st2 == 0 & (interim_prob_dlt_st2[c2 + 1] <= delta_Z2), 0, 1)
    sample_size_stopping_safety_st2 = ifelse(stopping_safety_ind_st2 == 0 & interim_prob_dlt_st2[c2 + 1] <= delta_Z2 & (sample_size_stopping_safety_st2 < (NN_1 + NN_2)), sample_size_stopping_safety_st2 + m2, sample_size_stopping_safety_st2)
    
    post_beta0  = median(ss$beta0)
    post_beta1  = median(ss$beta1)
    post_beta2  = median(ss$beta2)
    post_beta3  = median(ss$beta3)

    post_prob_eff = function(x){
      y_mtd = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3,theta,x)
      y_mtd = ifelse(y_mtd < 0, 0, ifelse(y_mtd > 1, 1, y_mtd))
      
      df_posterior = data.frame(beta0 = ss$beta0, beta1 = ss$beta1, beta2 = ss$beta2, 
                                beta3 = ss$beta3, x = x, y = y_mtd) %>%
        mutate(ind = 1:n()) %>%
        group_by(ind) %>%
        mutate(peff = peff(beta0,beta1,beta2,beta3,x,y),
               condition = ifelse(peff > p0, 1, 0)) %>%
        ungroup() %>%
        summarise(out = mean(condition))
      return(df_posterior$out)
    }
    post_prob_efficacy_estimated_MTD = as.numeric(lapply(x_grid, post_prob_eff))
    interim_prob_futility[c2] = max(post_prob_efficacy_estimated_MTD)
    
    #We assess whether we will enroll a new cohort or not based on futility stopping rule
    stopping_futility_ind = ifelse(stopping_futility_ind == 0 & (interim_prob_futility[c2] >= delta_E0), 0, 1)
    sample_size_stopping_futility = ifelse(stopping_futility_ind == 0 & interim_prob_futility[c2] >= delta_E0 & (sample_size_stopping_futility < (NN_1 + NN_2)), sample_size_stopping_futility + m2, sample_size_stopping_futility)
    
  }
  
  
  # ------------ Produce operating characteristics ------------- #
  
  #Proportion of patients in stage II allocated in doses with true prob of efficacy > theta
  
  proportion_allocation_efficacious_doses_df = data.frame(x = doseX2[(n2+1):NN_2],
                                                          y = doseY2[(n2+1):NN_2]) %>%
    mutate(peff = peff(true_beta0_st2, true_beta1_st2, true_beta2_st2, true_beta3_st2, x, y),
           indicator = ifelse(peff > p0, 1, 0)) %>%
    summarise(mean_allocation = mean(indicator))
  
  
  post_prob_efficacious_estimated_MTD = data.frame(x = x_grid) %>%
    mutate(y = twodimmtd2(post_rho00, post_rho01, post_rho10, post_alpha3, theta, x),
           post_prob_efficacious = post_prob_efficacy_estimated_MTD) %>%
    filter(y >= 0 & y <= 1) %>%
    mutate(obs_ind = 1:n()) %>%
    group_by(obs_ind) %>%
    mutate(true_peff = peff(true_beta0_st2, true_beta1_st2, true_beta2_st2, true_beta3_st2, x, y)) %>%
    ungroup()
  
  optimal_dose_combination = post_prob_efficacious_estimated_MTD %>%
    filter(post_prob_efficacious == max(post_prob_efficacious)) %>%
    dplyr::select(x,y,true_peff)
  
  #-----------------------------------------------------------------------------------------------#
  
  #Global output
  
  dt = list(Z.stg1 = Z,
            Z.stg2 = Z.stg2,
            Z.overall = c(Z,Z.stg2),
            NEff.stg1 = NEff.stg1,
            NEff.stg2 = NEff.stg2,
            doseX1 = doseX1,
            doseY1 = doseY1,
            doseX2 = doseX2[(n2+1):NN_2],
            doseY2 = doseY2[(n2+1):NN_2],
            proportion_allocation_efficacious_doses = proportion_allocation_efficacious_doses_df$mean_allocation,
            post_rho00 = post_rho00,
            post_rho10 = post_rho10,
            post_rho01 = post_rho01,
            post_alpha3 = post_alpha3,
            post_beta0_st2 = post_beta0_st2,
            post_beta1_st2 = post_beta1_st2,
            post_beta2_st2 = post_beta2_st2,
            post_beta3_st2 = post_beta3_st2,
            MTD_curve = twodimmtd2(true_rho00,true_rho01,true_rho10,true_alpha3,theta, seq(0,1,length.out = 101)),
            estimated_MTD_curve = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3,theta, seq(0,1,length.out = 101)),
            max_post_prob_efficacy_higher_theta_E = max(post_prob_efficacy_estimated_MTD),
            absolute_difference_rho00 = post_rho00 - true_rho00,
            absolute_difference_rho10 = post_rho10 - true_rho10,
            absolute_difference_rho01 = post_rho01 - true_rho01,
            absolute_difference_alpha3 = post_alpha3 - true_alpha3,
            absolute_difference_beta0_st2 = post_beta0_st2 - true_beta0_st2,
            absolute_difference_beta1_st2 = post_beta1_st2 - true_beta1_st2,
            absolute_difference_beta2_st2 = post_beta2_st2 - true_beta2_st2,
            absolute_difference_beta3_st2 = post_beta3_st2 - true_beta3_st2,
            interim_prob_dlt_st1 = interim_prob_dlt_st1,
            interim_prob_dlt_st2 = interim_prob_dlt_st2,
            interim_prob_futility = interim_prob_futility,
            sample_size_stopping_futility = sample_size_stopping_futility,
            sample_size_stopping_safety_st1 = sample_size_stopping_safety_st1,
            sample_size_stopping_safety_st2 = sample_size_stopping_safety_st2,
            optimal_dose_combination = optimal_dose_combination)
  
  dt
  
}








