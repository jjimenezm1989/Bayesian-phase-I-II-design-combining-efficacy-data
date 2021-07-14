#Add path to JAGS files.
setwd("")

library(rjags)
library(tidyverse)

rm(list=ls())

set.seed(1234)

#######################
#User-defined function#
#######################

# This function simulates a binary DLT responses for stages 1 and 2.
resp_safety_st1 = function(rho00 = trho00_st1, rho01 = trho01_st1, rho10 = trho10_st1, alpha3 = talpha3_st1, x, y){
  alpha0 = log(rho00/(1-rho00))
  alpha1 = log(rho10/(1-rho10)) - log(rho00/(1-rho00))
  alpha2 = log(rho01/(1-rho01)) - log(rho00/(1-rho00))
  p = 1/(1+exp(-alpha0-alpha1*x-alpha2*y-alpha3*x*y))
  u = runif(1)
  z = ifelse(u < p, 1, 0)
  return(z)
}
resp_safety_st2 = function(rho00 = trho00_st2, rho01 = trho01_st2, rho10 = trho10_st2, alpha3 = talpha3_st2, x, y){
  alpha0 = log(rho00/(1-rho00))
  alpha1 = log(rho10/(1-rho10)) - log(rho00/(1-rho00))
  alpha2 = log(rho01/(1-rho01)) - log(rho00/(1-rho00))
  p = 1/(1+exp(-alpha0-alpha1*x-alpha2*y-alpha3*x*y))
  u = runif(1)
  z = ifelse(u < p, 1, 0)
  return(z)
}

# This function simulates a binary efficacy responses.
resp_efficacy_st1 = function(beta0 = tbeta0_st1, beta1 = tbeta1_st1, beta2 = tbeta2_st1, beta3 = tbeta3_st1, x, y){
  p = exp(beta0+exp(beta1)*x+exp(beta2)*y+beta3*x*y) / (1 + exp(beta0+exp(beta1)*x+exp(beta2)*y+beta3*x*y))
  u = runif(1)
  z = ifelse(u < p, 1, 0)
  return(z)
}
resp_efficacy_st2 = function(beta0 = tbeta0_st2, beta1 = tbeta1_st2, beta2 = tbeta2_st2, beta3 = tbeta3_st2, x, y){
  p = exp(beta0+exp(beta1)*x+exp(beta2)*y+beta3*x*y) / (1 + exp(beta0+exp(beta1)*x+exp(beta2)*y+beta3*x*y))
  u = runif(1)
  z = ifelse(u < p, 1, 0)
  return(z)
}

# This function computes the probability of DLT at dose combination (x,y)
pdlt = function(rho00,rho01,rho10,alpha3,x,y){
  alpha0 = log(rho00/(1-rho00))
  alpha1 = log(rho10/(1-rho10)) - log(rho00/(1-rho00))
  alpha2 = log(rho01/(1-rho01)) - log(rho00/(1-rho00))
  p = 1/(1+exp(-alpha0-alpha1*x-alpha2*y-alpha3*x*y))
  return(p)
}

# This function computes the probability of efficacy at dose combination (x,y)
peff = function(beta0,beta1,beta2,beta3,x,y){
  p = exp(beta0+exp(beta1)*x+exp(beta2)*y+beta3*x*y) / (1 + exp(beta0+exp(beta1)*x+exp(beta2)*y+beta3*x*y))
  return(p)
}

# This function finds the set of MTD set
twodimmtd2 = function(rho00,rho01,rho10,alpha3,phi,x){
  
  #logistic(u) = exp(u)/(1+exp(u))
  #logistic(u)^-1 = logit(u) = log(u/(1-u))
  
  alpha0 = log(rho00/(1-rho00))
  alpha1 = log(rho10/(1-rho10)) - log(rho00/(1-rho00))
  alpha2 = log(rho01/(1-rho01)) - log(rho00/(1-rho00))
  
  a = (log(phi/(1-phi))-alpha0-alpha1*x)/(alpha2+alpha3*x)
  
  return(a)
}


#######################
#Main part of the code#
#######################

M_iter = 1000

# Controls size of jump
delta1 = 0.2
Xlow = 0
Ylow = 0

# Minimum and maximum doses available in the trial
# X dose is for Cabazitaxel
Xmin = 10  
Xmax = 25
# Y dose is for Cisplatin
Ymin = 50
Ymax = 100

lbx  =  (0 - Xmin)/(Xmax - Xmin)
lby  =  (0 - Ymin)/(Ymax - Ymin)

qq = 101
x_grid = seq(0,1,length.out = qq)

# Threshold for posterior probability of toxicity at minimum dose for stopping trial
delta = 0.8

#True model parameters

# True parameters for safety in stage 1
trho00_st1<-1e-7
trho01_st1<-0.2
trho10_st1<-0.2
talpha3_st1<-10
# True parameters for safety in stage 2
trho00_st2 = trho00_st1
trho01_st2 = trho01_st1
trho10_st2 = trho10_st1
talpha3_st2 = talpha3_st1
#True parameters for efficacy is stage 1
tbeta0_st1 = -5
tbeta1_st1 = 0.75
tbeta2_st1 = 1.51
tbeta3_st1 = 0.5
#True parameters for efficacy is stage 2
tbeta0_st2 = tbeta0_st1
tbeta1_st2 = tbeta1_st1
tbeta2_st2 = tbeta2_st1
tbeta3_st2 = tbeta3_st1
# Target Probability of DLT 
phi = 1/3

#Efficacy from the standard of care
null_p = 0.15

# Number of patients in a trial
NN_1 = 30
NN_2_fixed = 10
NN_2 = 30 - NN_2_fixed
m_cohort = NN_2/5

# Feasibility bound (EWOC)
alpha = numeric()
alphaa = seq(0.4,0.5,by=0.05)
for (ind in 1:3){
  alpha[ind] = alphaa[ind]
}
for (ind in 4:(NN_1-1)){
  alpha[ind] = alphaa[3]
}

vrho00 = numeric()
vrho01 = numeric()
vrho10 = numeric()
valpha3 = numeric()

vbeta0_st2 = numeric()
vbeta1_st2 = numeric()
vbeta2_st2 = numeric()
vbeta3_st2 = numeric()

vprob.ex = numeric()

#mcmc parameters

chains = 1
burn = 5000
mm = 2500

#Prior hyperparameters (stage 1)

a00 = 0.8
b00 = 7.2
a01 = 1.4
b01 = 5.6
a10 = 1.4
b10 = 5.6
a = 0.8
b = 0.0384

#How to define wMix
#100% non-exchangeability
#wMix = c(0,1)
#10% exchangeability
#wMix = c(1,0)
#50% exchangeability - 50% non-exchangeability
#wMix = c(0.5,0.5)

#Prior exchangeability and non-exchangeability belief for the design
w = 0
wMix = c(w,(1-w))

# Declaration of intermediate statistics
trcmtdx1 = numeric()
trcmtdx2 = numeric()
trcmtdy1 = numeric()
trcmtdy1 = numeric()

# Declaration of the output statistics
optimal_dose_x = numeric()
optimal_dose_y = numeric()
dosex_st1 = mat.or.vec(M_iter,NN_1)
dosey_st1 = mat.or.vec(M_iter,NN_1)
dosex_st2 = mat.or.vec(M_iter,(NN_2_fixed+NN_2))
dosey_st2 = mat.or.vec(M_iter,(NN_2_fixed+NN_2))
interim_prob_eff = mat.or.vec(M_iter,NN_2/m_cohort)
posterior_prob_efficacy = mat.or.vec(M_iter,qq)
dlt_st1 = mat.or.vec(M_iter,NN_1)
eff_st1 = mat.or.vec(M_iter,NN_1)
eff_st2 = mat.or.vec(M_iter,(NN_2_fixed+NN_2))
max_posterior_prob_efficacy = numeric()
target_dose = numeric()

for (kk in 1:M_iter) {
  
  cat("Trial simulation number:",kk,"\n")
  
  # data for first cohort of two patients- Z is DLT status
  
  # First two patients will be given Cisplatin 75 mg/m2 and Cabazitaxel 15 mg/m2
  doseX1 = c(0.3333333,0.3333333)
  doseY1 = c(0.5,0.5)
  
  # Simulate responses for patients 1 and 2 from logistic model
  
  zz1 = resp_safety_st1(x = doseX1[1], y = doseY1[1])
  zz2 = resp_safety_st1(x = doseX1[2], y = doseY1[2])
  ee1 = resp_efficacy_st1(x = doseX1[1], y = doseY1[1])
  ee2 = resp_efficacy_st1(x = doseX1[2], y = doseY1[2])
  Z = c(zz1,zz2)
  NEff.stg1 = c(ee1, ee2)
  Ncomb1 = 2 
  
  data_jags = list('Z'= Z, 'X'= doseX1, 'Y'= doseY1, 'phi'= phi, 'N'= Ncomb1, 'a00' = a00, 'b00' = b00, 
                   'a01' = a01, 'b01' = b01, 'a10' = a10, 'b10' = b10, 'a' = a, 'b' = b)
  
  # Get dose Y for patient 3 when dose X equals 0 and dose X for patient 4 when dose Y equals 0
  j=jags.model('stage1.bug.txt',data=data_jags,n.chains=chains,n.adapt=burn,quiet = TRUE)
  s=coda.samples(j,c('mtdx1','mtdy1','mtdx2','mtdy2'),mm,progress.bar="none")
  ss=as.data.frame(s[[1]])
  
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
  
  # Simulate responses for patients 3 and 4 from logistic model
  zz1 = resp_safety_st1(x = doseX1[3], y = doseY1[3])
  zz2 = resp_safety_st1(x = doseX1[4], y = doseY1[4])
  ee1 = resp_efficacy_st1(x = doseX1[3], y = doseY1[3])
  ee2 = resp_efficacy_st1(x = doseX1[4], y = doseY1[4])
  Z = c(Z,zz1,zz2)
  NEff.stg1 = c(NEff.stg1, ee1, ee2)
  Ncomb1 = 4
  
  data_jags = list('Z'= Z, 'X'= doseX1, 'Y'= doseY1, 'phi'= phi, 'N'= Ncomb1, 'a00' = a00, 'b00' = b00, 
                   'a01' = a01, 'b01' = b01, 'a10' = a10, 'b10' = b10, 'a' = a, 'b' = b)
  
  for (i in 3:(NN_1/2)) {
    
    j=jags.model('stage1.bug.txt',data=data_jags,n.chains=chains,n.adapt=burn,quiet = TRUE)
    s=coda.samples(j,c('rho00','rho01','rho10','alpha3','mtdx1','mtdy1','mtdx2','mtdy2'),mm,progress.bar="none")
    ss=as.data.frame(s[[1]])
    
    # #Calculating posterior probability of toxicity (UNCOMMENT LATER ON)
    # ss2 = ss %>%
    #   mutate(temp1 = pdlt(rho00,rho01,rho10,alpha3, 0, 0),
    #          temp2 = ifelse(temp1 > 0.43333, 1, 0))
    # 
    # postlow[kk,i] = as.numeric(ss2 %>% summarise(m = mean(temp2)))
    
    
    trcmtdx1 = ss$mtdx1[ss$mtdx1 > lbx]
    xx1 = max(0,quantile(trcmtdx1,alpha[i-1]))
    xx1 = min(xx1,1)
    if ((xx1 - doseX1[2*i-3]) > delta1)
      xx1 = doseX1[2*i-3]+delta1
    
    trcmtdx2 = ss$mtdx2[ss$mtdx2 > lbx]
    xx2 = max(0,quantile(trcmtdx2,alpha[i-1]))
    xx2 = min(xx2,1)
    if ((xx2 - doseX1[2*i-2]) > delta1)
      xx2 = doseX1[2*i-2]+delta1
    
    trcmtdy1 = ss$mtdy1[ss$mtdy1 > lby]
    yy1 = max(0,quantile(trcmtdy1,alpha[i-1]))
    yy1 = min(yy1,1)
    if ((yy1 - doseY1[2*i-3]) > delta1)
      yy1 = doseY1[2*i-3]+delta1
    
    trcmtdy2 = ss$mtdy2[ss$mtdy2 > lby]
    yy2 = max(0,quantile(trcmtdy2,alpha[i-1]))
    yy2 = min(yy2,1)
    if ((yy2 - doseY1[2*i-2]) > delta1)
      yy2 = doseY1[2*i-2]+delta1
    
    
    if (doseX1[2*i-3] == doseX1[2*i-5]) {
      
      doseX1 = c(doseX1,xx1,doseX1[2*i-2])
      doseY1 = c(doseY1,doseY1[2*i-3],yy2)} else {
        doseX1 = c(doseX1,doseX1[2*i-3],xx2)
        doseY1 = c(doseY1,yy1,doseY1[2*i-2])}
    
    Ncomb1 = Ncomb1+2
    
    zz1 = resp_safety_st1(x = doseX1[2*i-1], y = doseY1[2*i-1])
    zz2 = resp_safety_st1(x = doseX1[2*i-2], y = doseY1[2*i-2])
    ee1 = resp_efficacy_st1(x = doseX1[2*i-1], y = doseY1[2*i-1])
    ee2 = resp_efficacy_st1(x = doseX1[2*i-2], y = doseY1[2*i-2])
    Z = c(Z,zz1,zz2)
    NEff.stg1 = c(NEff.stg1, ee1, ee2)
  }
  
  post_rho00 = median(ss$rho00)
  post_rho01 = median(ss$rho01)
  post_rho10 = median(ss$rho10)
  post_alpha3 = median(ss$alpha3)
  
  
  #----------End of part 1
  
  #From now on, the MTD remains fixed since Stage 1 and Stage 2 populations are not the same. This is a 
  #limitation of this approach since most likely, the true MTD curve in Stage 2 may be slightly different.
  
  #We don't have any prior knowledge regarding efficacy on Stage 2. Therefore, we allocate 10 patients across the MTD
  #to collect efficacy data from Stage 2 on the MTD. After these 10 patients, we will do adaptive randomization.
  
  #We first see in which interval the MTD is defined for X.
  mtd_data_frame = data.frame(x = x_grid) %>%
    mutate(y = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3, phi = phi, x = x),
           y = ifelse(y < 0, 0, ifelse(y > 1, 1, y)))
  
  #We now know the in what interval X is defined. We allocate 10 patients equally spaced along the MTD
  xx1 = seq(0,1,length.out = NN_2_fixed)
  yy1 = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3, phi = phi, x = xx1)
  yy1 = ifelse(yy1 < 0, 0, ifelse(yy1 > 1, 1, yy1))
  
  #We now add the stage 2 data.
  doseX2 =  xx1
  doseY2 =  yy1
  Ncomb2 = NN_2_fixed
  NEff.stg2 = numeric()
  
  #We generate safety and efficacy outputs for the 10 patients placed along the MTD curve
  for(i in 1:NN_2_fixed){
    NEff.stg2[i] = resp_efficacy_st2(x = xx1[i], y = yy1[i])
  }
  
  data_jags = list('NEff.stg1' = NEff.stg1, 'doseX1' = doseX1, 'doseY1' = doseY1, 'Ncomb1' = Ncomb1,
                   'NEff.stg2' = NEff.stg2, 'doseX2' = doseX2, 'doseY2' = doseY2, 'Ncomb2' = Ncomb2,
                   'wMix' = wMix)
  
  j=jags.model('stage2_EXNEX.bug.txt',data= data_jags, n.chains=chains,n.adapt=burn,quiet = TRUE)
  s=coda.samples(j,c('beta02','beta12','beta22','beta32','prob.ex'),mm,progress.bar="none")
  ss=as.data.frame(s[[1]])
  
  post_beta0_st2 = median(ss$beta02)
  post_beta1_st2 = median(ss$beta12)
  post_beta2_st2 = median(ss$beta22)
  post_beta3_st2 = median(ss$beta32)
  post_prob.ex = mean(ss$`prob.ex[1]`)
  
  for(stage2 in 1:(NN_2/m_cohort)){
    
    #The MTD is updated after each cohort of m_cohort patients
    integrate_peff = function(x){
      y_mtd = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3,phi,x)
      y_mtd = ifelse(y_mtd < 0, 0, ifelse(y_mtd > 1, 1, y_mtd))
      p = peff(post_beta0_st2,post_beta1_st2,post_beta2_st2,post_beta3_st2,x,y_mtd)
      return(p)
      
    }
    
    ncst<-integrate(integrate_peff,0,1)$value
    
    #Get the maximum value of the density integrate_peff
    
    tempmax = dplyr::pull(data.frame(x_grid = x_grid) %>% 
                            mutate(y = integrate_peff(x_grid)) %>%
                            summarise(tempmax = max(y)))
    
    M<-(1/ncst)*tempmax
    
    #Sample from the standardized density based on the rejection-sampling principle.
    
    T = data.frame(sample.x = runif(600,0,1),
                   U = runif(600,0,1)) %>%
      mutate(rejection_boundary = (1/ncst)*integrate_peff(sample.x),
             density_sample_value = dunif(sample.x, 0, 1)*M*U,
             acceptance = ifelse(density_sample_value <= rejection_boundary, "Yes", "No")) %>%
      filter(acceptance == "Yes")
    
    T = T %>% dplyr::select(sample.x) %>% slice(1:m_cohort)
    
    T = T %>% mutate(y = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3,phi,sample.x),
                     y = ifelse(y < 0, 0, ifelse(y > 1, 1, y)))
    
    
    xx1 = T$sample.x
    yy1 = T$y
    
    for (ii in 1:m_cohort){
      
      ee1 = resp_efficacy_st2(tbeta0_st2, tbeta1_st2, tbeta2_st2, tbeta3_st2, x = xx1[ii], y = yy1[ii])
      
      NEff.stg2 = c(NEff.stg2,ee1)
      doseX2 =  c(doseX2,xx1[ii])
      doseY2 =  c(doseY2,yy1[ii])
      
    }
    
    Ncomb2 = length(doseX2)
    
    list('NEff.stg1' = NEff.stg1, 'doseX1' = doseX1, 'doseY1' = doseY1, 'Ncomb1' = Ncomb1,
         'NEff.stg2' = NEff.stg2, 'doseX2' = doseX2, 'doseY2' = doseY2, 'Ncomb2' = Ncomb2,
         'wMix' = wMix)
    
    j=jags.model('stage2_EXNEX.bug.txt',data= data_jags, n.chains=chains,n.adapt=burn,quiet = TRUE)
    s=coda.samples(j,c('beta02','beta12','beta22','beta32','prob.ex'),mm,progress.bar="none")
    ss=as.data.frame(s[[1]])
    
    post_beta0_st2 = median(ss$beta02)
    post_beta1_st2 = median(ss$beta12)
    post_beta2_st2 = median(ss$beta22)
    post_beta3_st2 = median(ss$beta32)
    post_prob.ex = mean(ss$`prob.ex[1]`)
    
    post_prob_eff = function(x){
      y_mtd = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3,phi,x)
      y_mtd = ifelse(y_mtd < 0, 0, ifelse(y_mtd > 1, 1, y_mtd))
      
      df_posterior = data.frame(beta02=ss$beta02,beta12=ss$beta12,beta22=ss$beta22,beta32=ss$beta32, x = x, y = y_mtd) %>%
        mutate(peff = peff(beta02,beta12,beta22,beta32,x,y),
               condition = ifelse(peff > null_p,1,0)) %>%
        summarise(out = mean(condition))
      return(df_posterior$out)
    }
    interim_prob_eff[kk,stage2] = max(as.numeric(lapply(x_grid, post_prob_eff)))
    
  }
  
  #Save trial doses recommended
  dosex_st1[kk,] = doseX1
  dosey_st1[kk,] = doseY1
  dosex_st2[kk,] = doseX2
  dosey_st2[kk,] = doseY2
  
  #Save trial safety and efficacy outcomes
  
  dlt_st1[kk,] = Z
  eff_st1[kk,] = NEff.stg1
  eff_st2[kk,] = NEff.stg2
  
  #Save trial posterior parameter estimates
  
  vrho00[kk] = post_rho00
  vrho01[kk] = post_rho01
  vrho10[kk] = post_rho10
  valpha3[kk] = post_alpha3
  vbeta0_st2[kk] = post_beta0_st2
  vbeta1_st2[kk] = post_beta1_st2
  vbeta2_st2[kk] = post_beta2_st2
  vbeta3_st2[kk] = post_beta3_st2
  vprob.ex[kk] = post_prob.ex
  
  post_prob_eff = function(x){
    y_mtd = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3,phi,x)
    y_mtd = ifelse(y_mtd < 0, 0, ifelse(y_mtd > 1, 1, y_mtd))
    
    df_posterior = data.frame(beta02=ss$beta02,beta12=ss$beta12,beta22=ss$beta22,beta32=ss$beta32, x = x, y = y_mtd) %>%
      mutate(peff = peff(beta02,beta12,beta22,beta32,x,y),
             condition = ifelse(peff > null_p,1,0)) %>%
      summarise(out = mean(condition))
    return(df_posterior$out)
  }
  posterior_prob_efficacy[kk,] = as.numeric(lapply(x_grid, post_prob_eff))
  max_posterior_prob_efficacy[kk] = max(posterior_prob_efficacy[kk,])
  target_dose[kk]<-which.max(posterior_prob_efficacy[kk,])
  
}


#END OF MAIN LOOP

#Files to save:

power = length(max_posterior_prob_efficacy[max_posterior_prob_efficacy > 0.4])/M_iter
write.csv(power,"/home/jimenjox/Research/EXNEX/EXNEX binary/R code/DLT Sc1 - EFF Sc1/power_w0.csv")
write.csv(posterior_prob_efficacy,"/home/jimenjox/Research/EXNEX/EXNEX binary/R code/DLT Sc1 - EFF Sc1/posterior_prob_efficacy_w0.csv")
bias = c(mean(vbeta1_st2 - tbeta1_st2), mean(vbeta2_st2 - tbeta2_st2))
mse = c(mean((vbeta1_st2 - tbeta1_st2)^2),mean((vbeta2_st2 - tbeta2_st2)^2))
write.csv(bias,"/home/jimenjox/Research/EXNEX/EXNEX binary/R code/DLT Sc1 - EFF Sc1/bias_w0.csv")
write.csv(mse,"/home/jimenjox/Research/EXNEX/EXNEX binary/R code/DLT Sc1 - EFF Sc1/mse_w0.csv")
write.csv(dosex_st2,"/home/jimenjox/Research/EXNEX/EXNEX binary/R code/DLT Sc1 - EFF Sc1/dosex_st2_w0.csv")
write.csv(dosey_st2,"/home/jimenjox/Research/EXNEX/EXNEX binary/R code/DLT Sc1 - EFF Sc1/dosey_st2_w0.csv")
write.csv(interim_prob_eff,"/home/jimenjox/Research/EXNEX/EXNEX binary/R code/DLT Sc1 - EFF Sc1/interim_prob_eff_w0.csv")
