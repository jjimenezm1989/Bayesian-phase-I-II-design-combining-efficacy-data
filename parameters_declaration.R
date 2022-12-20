#'Design parameters

# Target Probability of DLT 
theta = 1/3 

#Probability of standard treatment of care
p0 = 0.15

# Number of patients in a trial 
NN_1 = 30
NN_2 = 30
m1 = 2
m2 = 5
n2 = 10

# Controls size of jump
delta1 = 0.2
Xlow = 0
Ylow = 0

#Early stopping thresholds
delta_Z1 = 0.5 #safety stopping rule stage I
delta_Z2 = 0.9 #safety stopping rule stage II
delta_E0 = 0.1 #futility stopping rule

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

# Feasibility bound (EWOC)
alpha = numeric()
alphaa = seq(0.4,0.5,by=0.05)
for (ind in 1:3){
  alpha[ind] = alphaa[ind]
}
for (ind in 4:(NN_1-1)){
  alpha[ind] = alphaa[3]
}


# mcmc parameters
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

#Dose-toxicity and dose-efficacy model parameters

#Dose-toxicity scenario 1

  true_rho00_DLT_sc1  = 1e-7
  true_rho01_DLT_sc1  = 0.2
  true_rho10_DLT_sc1  = 0.2
  true_alpha3_DLT_sc1 = 10

#Dose-toxicity scenario 2

  true_rho00_DLT_sc2  = 0.001
  true_rho01_DLT_sc2  = 0.05
  true_rho10_DLT_sc2  = 0.05
  true_alpha3_DLT_sc2 = 10

#Dose-toxicity scenario 1, dose-efficacy scenario 1, H1
  
true_beta0_DLT_sc1_EFF_sc1_st2_H1 = -5
true_beta1_DLT_sc1_EFF_sc1_st2_H1 = 0.75
true_beta2_DLT_sc1_EFF_sc1_st2_H1 = 1.51
true_beta3_DLT_sc1_EFF_sc1_st2_H1 = 0.5
  
true_beta0_DLT_sc1_EFF_sc1_CA_st1_H1 = true_beta0_DLT_sc1_EFF_sc1_st2_H1
true_beta1_DLT_sc1_EFF_sc1_CA_st1_H1 = true_beta1_DLT_sc1_EFF_sc1_st2_H1
true_beta2_DLT_sc1_EFF_sc1_CA_st1_H1 = true_beta2_DLT_sc1_EFF_sc1_st2_H1
true_beta3_DLT_sc1_EFF_sc1_CA_st1_H1 = true_beta3_DLT_sc1_EFF_sc1_st2_H1
    
true_beta0_DLT_sc1_EFF_sc1_PA_st1_H1 = -5
true_beta1_DLT_sc1_EFF_sc1_PA_st1_H1 = 0.35
true_beta2_DLT_sc1_EFF_sc1_PA_st1_H1 = 1.11
true_beta3_DLT_sc1_EFF_sc1_PA_st1_H1 = 0.5
  
true_beta0_DLT_sc1_EFF_sc1_CD_st1_H1 = -5
true_beta1_DLT_sc1_EFF_sc1_CD_st1_H1 = 1.31
true_beta2_DLT_sc1_EFF_sc1_CD_st1_H1 = 0.75
true_beta3_DLT_sc1_EFF_sc1_CD_st1_H1 = 0.5
  
#Dose-toxicity scenario 1, dose-efficacy scenario 2, H1
  
true_beta0_DLT_sc1_EFF_sc2_st2_H1 = -5
true_beta1_DLT_sc1_EFF_sc2_st2_H1 = 1.5035
true_beta2_DLT_sc1_EFF_sc2_st2_H1 = 1.1
true_beta3_DLT_sc1_EFF_sc2_st2_H1 = 0.5
  
true_beta0_DLT_sc1_EFF_sc2_CA_st1_H1 = true_beta0_DLT_sc1_EFF_sc2_st2_H1
true_beta1_DLT_sc1_EFF_sc2_CA_st1_H1 = true_beta1_DLT_sc1_EFF_sc2_st2_H1
true_beta2_DLT_sc1_EFF_sc2_CA_st1_H1 = true_beta2_DLT_sc1_EFF_sc2_st2_H1
true_beta3_DLT_sc1_EFF_sc2_CA_st1_H1 = true_beta3_DLT_sc1_EFF_sc2_st2_H1
  
true_beta0_DLT_sc1_EFF_sc2_PA_st1_H1 = -5
true_beta1_DLT_sc1_EFF_sc2_PA_st1_H1 = 1.5
true_beta2_DLT_sc1_EFF_sc2_PA_st1_H1 = 0.2
true_beta3_DLT_sc1_EFF_sc2_PA_st1_H1 = 0.5
  
true_beta0_DLT_sc1_EFF_sc2_CD_st1_H1 = -8
true_beta1_DLT_sc1_EFF_sc2_CD_st1_H1 = -10
true_beta2_DLT_sc1_EFF_sc2_CD_st1_H1 = -10
true_beta3_DLT_sc1_EFF_sc2_CD_st1_H1 = 0
  
  
#Dose-toxicity scenario 2, dose-efficacy scenario 1, H1
  
true_beta0_DLT_sc2_EFF_sc1_st2_H1 = -6
true_beta1_DLT_sc2_EFF_sc1_st2_H1 = 1.2
true_beta2_DLT_sc2_EFF_sc1_st2_H1 = 1.623
true_beta3_DLT_sc2_EFF_sc1_st2_H1 = 0

true_beta0_DLT_sc2_EFF_sc1_CA_st1_H1 = true_beta0_DLT_sc2_EFF_sc1_st2_H1
true_beta1_DLT_sc2_EFF_sc1_CA_st1_H1 = true_beta1_DLT_sc2_EFF_sc1_st2_H1
true_beta2_DLT_sc2_EFF_sc1_CA_st1_H1 = true_beta2_DLT_sc2_EFF_sc1_st2_H1
true_beta3_DLT_sc2_EFF_sc1_CA_st1_H1 = true_beta3_DLT_sc2_EFF_sc1_st2_H1

true_beta0_DLT_sc2_EFF_sc1_PA_st1_H1 = -6
true_beta1_DLT_sc2_EFF_sc1_PA_st1_H1 = 1.4
true_beta2_DLT_sc2_EFF_sc1_PA_st1_H1 = 1.6
true_beta3_DLT_sc2_EFF_sc1_PA_st1_H1 = 0
  
true_beta0_DLT_sc2_EFF_sc1_CD_st1_H1 = -6
true_beta1_DLT_sc2_EFF_sc1_CD_st1_H1 = 1.623
true_beta2_DLT_sc2_EFF_sc1_CD_st1_H1 = 1.2
true_beta3_DLT_sc2_EFF_sc1_CD_st1_H1 = 0
  
#Dose-toxicity scenario 2, dose-efficacy scenario 2, H1
  
true_beta0_DLT_sc2_EFF_sc2_st2_H1 = -4
true_beta1_DLT_sc2_EFF_sc2_st2_H1 = 1.025
true_beta2_DLT_sc2_EFF_sc2_st2_H1 = 0.7
true_beta3_DLT_sc2_EFF_sc2_st2_H1 = 3
  
true_beta0_DLT_sc2_EFF_sc2_CA_st1_H1 = true_beta0_DLT_sc2_EFF_sc2_st2_H1
true_beta1_DLT_sc2_EFF_sc2_CA_st1_H1 = true_beta1_DLT_sc2_EFF_sc2_st2_H1
true_beta2_DLT_sc2_EFF_sc2_CA_st1_H1 = true_beta2_DLT_sc2_EFF_sc2_st2_H1
true_beta3_DLT_sc2_EFF_sc2_CA_st1_H1 = true_beta3_DLT_sc2_EFF_sc2_st2_H1

true_beta0_DLT_sc2_EFF_sc2_PA_st1_H1 = -4.5
true_beta1_DLT_sc2_EFF_sc2_PA_st1_H1 = 1.025
true_beta2_DLT_sc2_EFF_sc2_PA_st1_H1 = 0.7
true_beta3_DLT_sc2_EFF_sc2_PA_st1_H1 = 3
  
true_beta0_DLT_sc2_EFF_sc2_CD_st1_H1 = -8
true_beta1_DLT_sc2_EFF_sc2_CD_st1_H1 = -5
true_beta2_DLT_sc2_EFF_sc2_CD_st1_H1 = -5
true_beta3_DLT_sc2_EFF_sc2_CD_st1_H1 = 27




#Dose-toxicity scenario 1, dose-efficacy scenario 1, H0

true_beta0_DLT_sc1_EFF_sc1_st2_H0 = -4
true_beta1_DLT_sc1_EFF_sc1_st2_H0 = -2
true_beta2_DLT_sc1_EFF_sc1_st2_H0 = 0.8
true_beta3_DLT_sc1_EFF_sc1_st2_H0 = 0.5

true_beta0_DLT_sc1_EFF_sc1_CA_st1_H0 = true_beta0_DLT_sc1_EFF_sc1_st2_H0
true_beta1_DLT_sc1_EFF_sc1_CA_st1_H0 = true_beta1_DLT_sc1_EFF_sc1_st2_H0
true_beta2_DLT_sc1_EFF_sc1_CA_st1_H0 = true_beta2_DLT_sc1_EFF_sc1_st2_H0
true_beta3_DLT_sc1_EFF_sc1_CA_st1_H0 = true_beta3_DLT_sc1_EFF_sc1_st2_H0

true_beta0_DLT_sc1_EFF_sc1_PA_st1_H0 = -4
true_beta1_DLT_sc1_EFF_sc1_PA_st1_H0 = -2
true_beta2_DLT_sc1_EFF_sc1_PA_st1_H0 = 0.1
true_beta3_DLT_sc1_EFF_sc1_PA_st1_H0 = 1

true_beta0_DLT_sc1_EFF_sc1_CD_st1_H0 = -4.5
true_beta1_DLT_sc1_EFF_sc1_CD_st1_H0 = 1.4
true_beta2_DLT_sc1_EFF_sc1_CD_st1_H0 = 1
true_beta3_DLT_sc1_EFF_sc1_CD_st1_H0 = 0.5


#Dose-toxicity scenario 1, dose-efficacy scenario 2, H0

true_beta0_DLT_sc1_EFF_sc2_st2_H0 = -6.36
true_beta1_DLT_sc1_EFF_sc2_st2_H0 = 1.5035
true_beta2_DLT_sc1_EFF_sc2_st2_H0 = 1.1
true_beta3_DLT_sc1_EFF_sc2_st2_H0 = 0.5

true_beta0_DLT_sc1_EFF_sc2_CA_st1_H0 = true_beta0_DLT_sc1_EFF_sc2_st2_H0
true_beta1_DLT_sc1_EFF_sc2_CA_st1_H0 = true_beta1_DLT_sc1_EFF_sc2_st2_H0
true_beta2_DLT_sc1_EFF_sc2_CA_st1_H0 = true_beta2_DLT_sc1_EFF_sc2_st2_H0
true_beta3_DLT_sc1_EFF_sc2_CA_st1_H0 = true_beta3_DLT_sc1_EFF_sc2_st2_H0

true_beta0_DLT_sc1_EFF_sc2_PA_st1_H0 = -6.36
true_beta1_DLT_sc1_EFF_sc2_PA_st1_H0 = 1.65
true_beta2_DLT_sc1_EFF_sc2_PA_st1_H0 = 1.1
true_beta3_DLT_sc1_EFF_sc2_PA_st1_H0 = 0.5

true_beta0_DLT_sc1_EFF_sc2_CD_st1_H0 = -1.5
true_beta1_DLT_sc1_EFF_sc2_CD_st1_H0 = -1
true_beta2_DLT_sc1_EFF_sc2_CD_st1_H0 = -1
true_beta3_DLT_sc1_EFF_sc2_CD_st1_H0 = 0.25


#Dose-toxicity scenario 2, dose-efficacy scenario 1, H0

true_beta0_DLT_sc2_EFF_sc1_st2_H0 = -6
true_beta1_DLT_sc2_EFF_sc1_st2_H0 = 1.1
true_beta2_DLT_sc2_EFF_sc1_st2_H0 = 1.323
true_beta3_DLT_sc2_EFF_sc1_st2_H0 = 0

true_beta0_DLT_sc2_EFF_sc1_CA_st1_H0 = true_beta0_DLT_sc2_EFF_sc1_st2_H0
true_beta1_DLT_sc2_EFF_sc1_CA_st1_H0 = true_beta1_DLT_sc2_EFF_sc1_st2_H0
true_beta2_DLT_sc2_EFF_sc1_CA_st1_H0 = true_beta2_DLT_sc2_EFF_sc1_st2_H0
true_beta3_DLT_sc2_EFF_sc1_CA_st1_H0 = true_beta3_DLT_sc2_EFF_sc1_st2_H0

true_beta0_DLT_sc2_EFF_sc1_PA_st1_H0 = -7
true_beta1_DLT_sc2_EFF_sc1_PA_st1_H0 = 1.1
true_beta2_DLT_sc2_EFF_sc1_PA_st1_H0 = 1.4
true_beta3_DLT_sc2_EFF_sc1_PA_st1_H0 = 0

true_beta0_DLT_sc2_EFF_sc1_CD_st1_H0 = -6
true_beta1_DLT_sc2_EFF_sc1_CD_st1_H0 = 1.622
true_beta2_DLT_sc2_EFF_sc1_CD_st1_H0 = 1.2
true_beta3_DLT_sc2_EFF_sc1_CD_st1_H0 = 0

#Dose-toxicity scenario 2, dose-efficacy scenario 2, H0

true_beta0_DLT_sc2_EFF_sc2_st2_H0 = -5.35
true_beta1_DLT_sc2_EFF_sc2_st2_H0 = 1.025
true_beta2_DLT_sc2_EFF_sc2_st2_H0 = 0.7
true_beta3_DLT_sc2_EFF_sc2_st2_H0 = 3

true_beta0_DLT_sc2_EFF_sc2_CA_st1_H0 = true_beta0_DLT_sc2_EFF_sc2_st2_H0
true_beta1_DLT_sc2_EFF_sc2_CA_st1_H0 = true_beta1_DLT_sc2_EFF_sc2_st2_H0
true_beta2_DLT_sc2_EFF_sc2_CA_st1_H0 = true_beta2_DLT_sc2_EFF_sc2_st2_H0
true_beta3_DLT_sc2_EFF_sc2_CA_st1_H0 = true_beta3_DLT_sc2_EFF_sc2_st2_H0

true_beta0_DLT_sc2_EFF_sc2_PA_st1_H0 = -5
true_beta1_DLT_sc2_EFF_sc2_PA_st1_H0 = 1.1
true_beta2_DLT_sc2_EFF_sc2_PA_st1_H0 = -1
true_beta3_DLT_sc2_EFF_sc2_PA_st1_H0 = 1

true_beta0_DLT_sc2_EFF_sc2_CD_st1_H0 = -3.54
true_beta1_DLT_sc2_EFF_sc2_CD_st1_H0 = 0.5
true_beta2_DLT_sc2_EFF_sc2_CD_st1_H0 = 1
true_beta3_DLT_sc2_EFF_sc2_CD_st1_H0 = 1
  
  