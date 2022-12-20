setwd("/home/jimenjox/Research/MAC phase I-II/R code/")

library(rjags)
library(tidyverse)
library(scales)
library(mvtnorm)
library(truncnorm)
library(data.table)
library(ggpubr)

rm(list=ls())

########################
#Scenario specification#
########################

#Seet seed based on date we run the simulations (07/02/2022)
time_seed = 07022022
set.seed(time_seed)

# Number of Trials to simulate
M_iter = 1000

#Choose the toxicity (DLT) scenario number and efficacy (EFF) scenario number
DLT_scenario_variables = "_DLT_sc2"
EFF_scenario_variables = "_EFF_sc1"
DLT_scenario_full_description = "Dose-toxicity scenario 2"
EFF_scenario_full_description = "Dose-efficacy scenario 1"
hypothesis = "_H0"
st1 = "_st1"
st2 = "_st2"
agreement = "_PA"
agreement_text = "PA"

#Choose between "Power" (if under H1) or "Type-I error (if under H0)".
power_or_typeIerror = "Type-I error"

#######################
#User-defined function#
#######################

source("generate_bivariate_outcome.R")
source("pdlt.R")
source("peff.R")
source("twodimmtd2.R")
source("parameters_declaration.R")
source("trial_simulation.R")

path_JAGS = getwd()

#######################
#Main part of the code#
#######################

true_rho00 = eval(parse(text = paste0("true_rho00",DLT_scenario_variables)))
true_rho01 = eval(parse(text = paste0("true_rho01",DLT_scenario_variables)))
true_rho10 = eval(parse(text = paste0("true_rho10",DLT_scenario_variables)))
true_alpha3 = eval(parse(text = paste0("true_alpha3",DLT_scenario_variables)))

true_beta0_st1 = eval(parse(text = paste0("true_beta0",DLT_scenario_variables,EFF_scenario_variables,agreement,st1,hypothesis)))
true_beta1_st1 = eval(parse(text = paste0("true_beta1",DLT_scenario_variables,EFF_scenario_variables,agreement,st1,hypothesis)))
true_beta2_st1 = eval(parse(text = paste0("true_beta2",DLT_scenario_variables,EFF_scenario_variables,agreement,st1,hypothesis)))
true_beta3_st1 = eval(parse(text = paste0("true_beta3",DLT_scenario_variables,EFF_scenario_variables,agreement,st1,hypothesis)))

true_beta0_st2 = eval(parse(text = paste0("true_beta0",DLT_scenario_variables,EFF_scenario_variables,st2,hypothesis)))
true_beta1_st2 = eval(parse(text = paste0("true_beta1",DLT_scenario_variables,EFF_scenario_variables,st2,hypothesis)))
true_beta2_st2 = eval(parse(text = paste0("true_beta2",DLT_scenario_variables,EFF_scenario_variables,st2,hypothesis)))
true_beta3_st2 = eval(parse(text = paste0("true_beta3",DLT_scenario_variables,EFF_scenario_variables,st2,hypothesis)))

start.time <- Sys.time()

omega = "_w0"
omega_num = 0
do.call("<-",list(paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                  replicate(M_iter, trial_simulation(true_rho00 = true_rho00,
                                                     true_rho01 = true_rho01,
                                                     true_rho10 = true_rho10,
                                                     true_alpha3 = true_alpha3,
                                                     true_beta0_st1 = true_beta0_st1,
                                                     true_beta1_st1 = true_beta1_st1,
                                                     true_beta2_st1 = true_beta2_st1,
                                                     true_beta3_st1 = true_beta3_st1,
                                                     true_beta0_st2 = true_beta0_st2,
                                                     true_beta1_st2 = true_beta1_st2,
                                                     true_beta2_st2 = true_beta2_st2,
                                                     true_beta3_st2 = true_beta3_st2,
                                                     path_JAGS = path_JAGS, 
                                                     w = omega_num), simplify = FALSE)))

omega = "_w025"
omega_num = 0.25
do.call("<-",list(paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                  replicate(M_iter, trial_simulation(true_rho00 = true_rho00,
                                                     true_rho01 = true_rho01,
                                                     true_rho10 = true_rho10,
                                                     true_alpha3 = true_alpha3,
                                                     true_beta0_st1 = true_beta0_st1,
                                                     true_beta1_st1 = true_beta1_st1,
                                                     true_beta2_st1 = true_beta2_st1,
                                                     true_beta3_st1 = true_beta3_st1,
                                                     true_beta0_st2 = true_beta0_st2,
                                                     true_beta1_st2 = true_beta1_st2,
                                                     true_beta2_st2 = true_beta2_st2,
                                                     true_beta3_st2 = true_beta3_st2,
                                                     path_JAGS = path_JAGS, 
                                                     w = omega_num), simplify = FALSE)))

omega = "_w050"
omega_num = 0.5
do.call("<-",list(paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                  replicate(M_iter, trial_simulation(true_rho00 = true_rho00,
                                                     true_rho01 = true_rho01,
                                                     true_rho10 = true_rho10,
                                                     true_alpha3 = true_alpha3,
                                                     true_beta0_st1 = true_beta0_st1,
                                                     true_beta1_st1 = true_beta1_st1,
                                                     true_beta2_st1 = true_beta2_st1,
                                                     true_beta3_st1 = true_beta3_st1,
                                                     true_beta0_st2 = true_beta0_st2,
                                                     true_beta1_st2 = true_beta1_st2,
                                                     true_beta2_st2 = true_beta2_st2,
                                                     true_beta3_st2 = true_beta3_st2,
                                                     path_JAGS = path_JAGS, 
                                                     w = omega_num), simplify = FALSE)))

omega = "_w075"
omega_num = 0.75
do.call("<-",list(paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                  replicate(M_iter, trial_simulation(true_rho00 = true_rho00,
                                                     true_rho01 = true_rho01,
                                                     true_rho10 = true_rho10,
                                                     true_alpha3 = true_alpha3,
                                                     true_beta0_st1 = true_beta0_st1,
                                                     true_beta1_st1 = true_beta1_st1,
                                                     true_beta2_st1 = true_beta2_st1,
                                                     true_beta3_st1 = true_beta3_st1,
                                                     true_beta0_st2 = true_beta0_st2,
                                                     true_beta1_st2 = true_beta1_st2,
                                                     true_beta2_st2 = true_beta2_st2,
                                                     true_beta3_st2 = true_beta3_st2,
                                                     path_JAGS = path_JAGS, 
                                                     w = omega_num), simplify = FALSE)))


omega = "_w1"
omega_num = 1
do.call("<-",list(paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                  replicate(M_iter, trial_simulation(true_rho00 = true_rho00,
                                                     true_rho01 = true_rho01,
                                                     true_rho10 = true_rho10,
                                                     true_alpha3 = true_alpha3,
                                                     true_beta0_st1 = true_beta0_st1,
                                                     true_beta1_st1 = true_beta1_st1,
                                                     true_beta2_st1 = true_beta2_st1,
                                                     true_beta3_st1 = true_beta3_st1,
                                                     true_beta0_st2 = true_beta0_st2,
                                                     true_beta1_st2 = true_beta1_st2,
                                                     true_beta2_st2 = true_beta2_st2,
                                                     true_beta3_st2 = true_beta3_st2,
                                                     path_JAGS = path_JAGS, 
                                                     w = omega_num), simplify = FALSE)))

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

save.image(paste0(path_JAGS,"/environment_sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".RData"))

#--------DERIVE OPERATING CHARACTERISTICS AND SAVE THEM IN .CSV FILES#--------#

setwd("/home/jimenjox/Research/MAC phase I-II/R code/")

source("get_bias_MSE.R")
source("get_power.R")
source("get_posterior_probabilities_stopping_rules.R")
source("get_proportion_correct_patient_allocation.R")
source("get_dose_combinations.R")
source("get_MTD_curve.R")
source("get_average_DLT_rate.R")
source("get_sample_size_stopping_rules.R")
source("get_optimal_dose_combination.R")
source("get_proportion_optimal_dose_combination_above_efficacy_threshold.R")

setwd("/home/jimenjox/Research/MAC phase I-II/R code/intermediate_outputs/")



#OMEGA = 0

omega = "_w0"
omega_num = 0

#Bias and MSE model parameters
write.csv(do.call("<-",list(paste0("bias_MSE_parameters",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_bias_MSE(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                         DLT_scenario = DLT_scenario_full_description,
                                         EFF_scenario = EFF_scenario_full_description,
                                         omega = omega_num,
                                         agreement = agreement_text))),
          paste0("bias_MSE_parameters",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))

#optimal dose combination
write.csv(do.call("<-",list(paste0("optimal_dose_combination",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_optimal_dose_combination(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                         DLT_scenario = DLT_scenario_full_description,
                                                         EFF_scenario = EFF_scenario_full_description,
                                                         omega = omega_num,
                                                         agreement = agreement_text))),
          paste0("optimal_dose_combination",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))

#Power / type-I error

write.csv(do.call("<-",list(paste0("power_or_typeIerror",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_power(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                      type = power_or_typeIerror,
                                      DLT_scenario = DLT_scenario_full_description,
                                      EFF_scenario = EFF_scenario_full_description,
                                      omega = omega_num,
                                      agreement = agreement_text))),
          paste0("power_or_typeIerror",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))


#Stopping rules
write.csv(do.call("<-",list(paste0("stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_posterior_probabilities_stopping_rules(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                       delta_Z1 = delta_Z1,
                                                                       delta_Z2 = delta_Z2,
                                                                       delta_E0 = delta_E0,
                                                                       DLT_scenario = DLT_scenario_full_description,
                                                                       EFF_scenario = EFF_scenario_full_description,
                                                                       omega = omega_num,
                                                                       agreement = agreement_text))),
          paste0("stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#Sample size stopping rules
write.csv(do.call("<-",list(paste0("sample_size_stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_sample_size_stopping_rules(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                           DLT_scenario = DLT_scenario_full_description,
                                                           EFF_scenario = EFF_scenario_full_description,
                                                           omega = omega_num,
                                                           agreement = agreement_text))),
          paste0("sample_size_stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#Proportion patient optimal dose combinations with true P(pi_E) > theta_E
write.csv(do.call("<-",list(paste0("optimal_dose_combination_above_threshold",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_proportion_optimal_dose_combination_above_efficacy_threshold(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                                             DLT_scenario = DLT_scenario_full_description,
                                                                                             EFF_scenario = EFF_scenario_full_description,
                                                                                             omega = omega_num,
                                                                                             agreement = agreement_text))),
          paste0("optimal_dose_combination_above_threshold",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))




#Proportion patient allocation is dose combination with true P(pi_E) > theta_E
write.csv(do.call("<-",list(paste0("proportion_allocation",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_proportion_correct_patient_allocation(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                      DLT_scenario = DLT_scenario_full_description,
                                                                      EFF_scenario = EFF_scenario_full_description,
                                                                      omega = omega_num,
                                                                      agreement = agreement_text))),
          paste0("proportion_allocation",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))


#Dose combinations to which patients were allocated
write.csv(do.call("<-",list(paste0("dose_combinations",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_dose_combinations(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                  DLT_scenario = DLT_scenario_full_description,
                                                  EFF_scenario = EFF_scenario_full_description,
                                                  omega = omega_num,
                                                  agreement = agreement_text))),
          paste0("dose_combinations",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))




# MTD curve
write.csv(do.call("<-",list(paste0("MTD_curve",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_MTD_curve(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                          DLT_scenario = DLT_scenario_full_description,
                                          EFF_scenario = EFF_scenario_full_description,
                                          omega = omega_num,
                                          agreement = agreement_text))),
          paste0("MTD_curve",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#DLT rate and proportion of trials with DLT rate above theta_Z + 0.1
write.csv(do.call("<-",list(paste0("DLT_rate",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_average_DLT_rate(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                 DLT_scenario = DLT_scenario_full_description,
                                                 EFF_scenario = EFF_scenario_full_description,
                                                 omega = omega_num,
                                                 agreement = agreement_text))),
          paste0("DLT_rate",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))


#OMEGA = 0.25

omega = "_w025"
omega_num = 0.25

#Bias and MSE model parameters
write.csv(do.call("<-",list(paste0("bias_MSE_parameters",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_bias_MSE(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                         DLT_scenario = DLT_scenario_full_description,
                                         EFF_scenario = EFF_scenario_full_description,
                                         omega = omega_num,
                                         agreement = agreement_text))),
          paste0("bias_MSE_parameters",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))

#optimal dose combination
write.csv(do.call("<-",list(paste0("optimal_dose_combination",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_optimal_dose_combination(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                         DLT_scenario = DLT_scenario_full_description,
                                                         EFF_scenario = EFF_scenario_full_description,
                                                         omega = omega_num,
                                                         agreement = agreement_text))),
          paste0("optimal_dose_combination",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))

#Power / type-I error

write.csv(do.call("<-",list(paste0("power_or_typeIerror",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_power(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                      type = power_or_typeIerror,
                                      DLT_scenario = DLT_scenario_full_description,
                                      EFF_scenario = EFF_scenario_full_description,
                                      omega = omega_num,
                                      agreement = agreement_text))),
          paste0("power_or_typeIerror",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))


#Stopping rules
write.csv(do.call("<-",list(paste0("stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_posterior_probabilities_stopping_rules(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                       delta_Z1 = delta_Z1,
                                                                       delta_Z2 = delta_Z2,
                                                                       delta_E0 = delta_E0,
                                                                       DLT_scenario = DLT_scenario_full_description,
                                                                       EFF_scenario = EFF_scenario_full_description,
                                                                       omega = omega_num,
                                                                       agreement = agreement_text))),
          paste0("stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#Sample size stopping rules
write.csv(do.call("<-",list(paste0("sample_size_stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_sample_size_stopping_rules(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                           DLT_scenario = DLT_scenario_full_description,
                                                           EFF_scenario = EFF_scenario_full_description,
                                                           omega = omega_num,
                                                           agreement = agreement_text))),
          paste0("sample_size_stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#Proportion patient optimal dose combinations with true P(pi_E) > theta_E
write.csv(do.call("<-",list(paste0("optimal_dose_combination_above_threshold",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_proportion_optimal_dose_combination_above_efficacy_threshold(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                                             DLT_scenario = DLT_scenario_full_description,
                                                                                             EFF_scenario = EFF_scenario_full_description,
                                                                                             omega = omega_num,
                                                                                             agreement = agreement_text))),
          paste0("optimal_dose_combination_above_threshold",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))




#Proportion patient allocation is dose combination with true P(pi_E) > theta_E
write.csv(do.call("<-",list(paste0("proportion_allocation",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_proportion_correct_patient_allocation(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                      DLT_scenario = DLT_scenario_full_description,
                                                                      EFF_scenario = EFF_scenario_full_description,
                                                                      omega = omega_num,
                                                                      agreement = agreement_text))),
          paste0("proportion_allocation",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))


#Dose combinations to which patients were allocated
write.csv(do.call("<-",list(paste0("dose_combinations",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_dose_combinations(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                  DLT_scenario = DLT_scenario_full_description,
                                                  EFF_scenario = EFF_scenario_full_description,
                                                  omega = omega_num,
                                                  agreement = agreement_text))),
          paste0("dose_combinations",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))




# MTD curve
write.csv(do.call("<-",list(paste0("MTD_curve",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_MTD_curve(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                          DLT_scenario = DLT_scenario_full_description,
                                          EFF_scenario = EFF_scenario_full_description,
                                          omega = omega_num,
                                          agreement = agreement_text))),
          paste0("MTD_curve",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#DLT rate and proportion of trials with DLT rate above theta_Z + 0.1
write.csv(do.call("<-",list(paste0("DLT_rate",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_average_DLT_rate(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                 DLT_scenario = DLT_scenario_full_description,
                                                 EFF_scenario = EFF_scenario_full_description,
                                                 omega = omega_num,
                                                 agreement = agreement_text))),
          paste0("DLT_rate",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))


#OMEGA = 0.5

omega = "_w050"
omega_num = 0.5

#Bias and MSE model parameters
write.csv(do.call("<-",list(paste0("bias_MSE_parameters",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_bias_MSE(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                         DLT_scenario = DLT_scenario_full_description,
                                         EFF_scenario = EFF_scenario_full_description,
                                         omega = omega_num,
                                         agreement = agreement_text))),
          paste0("bias_MSE_parameters",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))

#optimal dose combination
write.csv(do.call("<-",list(paste0("optimal_dose_combination",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_optimal_dose_combination(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                         DLT_scenario = DLT_scenario_full_description,
                                                         EFF_scenario = EFF_scenario_full_description,
                                                         omega = omega_num,
                                                         agreement = agreement_text))),
          paste0("optimal_dose_combination",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))

#Power / type-I error

write.csv(do.call("<-",list(paste0("power_or_typeIerror",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_power(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                      type = power_or_typeIerror,
                                      DLT_scenario = DLT_scenario_full_description,
                                      EFF_scenario = EFF_scenario_full_description,
                                      omega = omega_num,
                                      agreement = agreement_text))),
          paste0("power_or_typeIerror",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))


#Stopping rules
write.csv(do.call("<-",list(paste0("stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_posterior_probabilities_stopping_rules(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                       delta_Z1 = delta_Z1,
                                                                       delta_Z2 = delta_Z2,
                                                                       delta_E0 = delta_E0,
                                                                       DLT_scenario = DLT_scenario_full_description,
                                                                       EFF_scenario = EFF_scenario_full_description,
                                                                       omega = omega_num,
                                                                       agreement = agreement_text))),
          paste0("stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#Sample size stopping rules
write.csv(do.call("<-",list(paste0("sample_size_stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_sample_size_stopping_rules(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                           DLT_scenario = DLT_scenario_full_description,
                                                           EFF_scenario = EFF_scenario_full_description,
                                                           omega = omega_num,
                                                           agreement = agreement_text))),
          paste0("sample_size_stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#Proportion patient optimal dose combinations with true P(pi_E) > theta_E
write.csv(do.call("<-",list(paste0("optimal_dose_combination_above_threshold",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_proportion_optimal_dose_combination_above_efficacy_threshold(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                                             DLT_scenario = DLT_scenario_full_description,
                                                                                             EFF_scenario = EFF_scenario_full_description,
                                                                                             omega = omega_num,
                                                                                             agreement = agreement_text))),
          paste0("optimal_dose_combination_above_threshold",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))




#Proportion patient allocation is dose combination with true P(pi_E) > theta_E
write.csv(do.call("<-",list(paste0("proportion_allocation",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_proportion_correct_patient_allocation(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                      DLT_scenario = DLT_scenario_full_description,
                                                                      EFF_scenario = EFF_scenario_full_description,
                                                                      omega = omega_num,
                                                                      agreement = agreement_text))),
          paste0("proportion_allocation",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))


#Dose combinations to which patients were allocated
write.csv(do.call("<-",list(paste0("dose_combinations",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_dose_combinations(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                  DLT_scenario = DLT_scenario_full_description,
                                                  EFF_scenario = EFF_scenario_full_description,
                                                  omega = omega_num,
                                                  agreement = agreement_text))),
          paste0("dose_combinations",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))




# MTD curve
write.csv(do.call("<-",list(paste0("MTD_curve",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_MTD_curve(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                          DLT_scenario = DLT_scenario_full_description,
                                          EFF_scenario = EFF_scenario_full_description,
                                          omega = omega_num,
                                          agreement = agreement_text))),
          paste0("MTD_curve",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#DLT rate and proportion of trials with DLT rate above theta_Z + 0.1
write.csv(do.call("<-",list(paste0("DLT_rate",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_average_DLT_rate(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                 DLT_scenario = DLT_scenario_full_description,
                                                 EFF_scenario = EFF_scenario_full_description,
                                                 omega = omega_num,
                                                 agreement = agreement_text))),
          paste0("DLT_rate",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))


#OMEGA = 0.75

omega = "_w075"
omega_num = 0.75

#Bias and MSE model parameters
write.csv(do.call("<-",list(paste0("bias_MSE_parameters",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_bias_MSE(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                         DLT_scenario = DLT_scenario_full_description,
                                         EFF_scenario = EFF_scenario_full_description,
                                         omega = omega_num,
                                         agreement = agreement_text))),
          paste0("bias_MSE_parameters",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))

#optimal dose combination
write.csv(do.call("<-",list(paste0("optimal_dose_combination",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_optimal_dose_combination(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                         DLT_scenario = DLT_scenario_full_description,
                                                         EFF_scenario = EFF_scenario_full_description,
                                                         omega = omega_num,
                                                         agreement = agreement_text))),
          paste0("optimal_dose_combination",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))

#Power / type-I error

write.csv(do.call("<-",list(paste0("power_or_typeIerror",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_power(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                      type = power_or_typeIerror,
                                      DLT_scenario = DLT_scenario_full_description,
                                      EFF_scenario = EFF_scenario_full_description,
                                      omega = omega_num,
                                      agreement = agreement_text))),
          paste0("power_or_typeIerror",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))


#Stopping rules
write.csv(do.call("<-",list(paste0("stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_posterior_probabilities_stopping_rules(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                       delta_Z1 = delta_Z1,
                                                                       delta_Z2 = delta_Z2,
                                                                       delta_E0 = delta_E0,
                                                                       DLT_scenario = DLT_scenario_full_description,
                                                                       EFF_scenario = EFF_scenario_full_description,
                                                                       omega = omega_num,
                                                                       agreement = agreement_text))),
          paste0("stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#Sample size stopping rules
write.csv(do.call("<-",list(paste0("sample_size_stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_sample_size_stopping_rules(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                           DLT_scenario = DLT_scenario_full_description,
                                                           EFF_scenario = EFF_scenario_full_description,
                                                           omega = omega_num,
                                                           agreement = agreement_text))),
          paste0("sample_size_stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#Proportion patient optimal dose combinations with true P(pi_E) > theta_E
write.csv(do.call("<-",list(paste0("optimal_dose_combination_above_threshold",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_proportion_optimal_dose_combination_above_efficacy_threshold(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                                             DLT_scenario = DLT_scenario_full_description,
                                                                                             EFF_scenario = EFF_scenario_full_description,
                                                                                             omega = omega_num,
                                                                                             agreement = agreement_text))),
          paste0("optimal_dose_combination_above_threshold",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))




#Proportion patient allocation is dose combination with true P(pi_E) > theta_E
write.csv(do.call("<-",list(paste0("proportion_allocation",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_proportion_correct_patient_allocation(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                      DLT_scenario = DLT_scenario_full_description,
                                                                      EFF_scenario = EFF_scenario_full_description,
                                                                      omega = omega_num,
                                                                      agreement = agreement_text))),
          paste0("proportion_allocation",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))


#Dose combinations to which patients were allocated
write.csv(do.call("<-",list(paste0("dose_combinations",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_dose_combinations(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                  DLT_scenario = DLT_scenario_full_description,
                                                  EFF_scenario = EFF_scenario_full_description,
                                                  omega = omega_num,
                                                  agreement = agreement_text))),
          paste0("dose_combinations",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))




# MTD curve
write.csv(do.call("<-",list(paste0("MTD_curve",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_MTD_curve(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                          DLT_scenario = DLT_scenario_full_description,
                                          EFF_scenario = EFF_scenario_full_description,
                                          omega = omega_num,
                                          agreement = agreement_text))),
          paste0("MTD_curve",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#DLT rate and proportion of trials with DLT rate above theta_Z + 0.1
write.csv(do.call("<-",list(paste0("DLT_rate",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_average_DLT_rate(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                 DLT_scenario = DLT_scenario_full_description,
                                                 EFF_scenario = EFF_scenario_full_description,
                                                 omega = omega_num,
                                                 agreement = agreement_text))),
          paste0("DLT_rate",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#OMEGA = 1

omega = "_w1"
omega_num = 1

#Bias and MSE model parameters
write.csv(do.call("<-",list(paste0("bias_MSE_parameters",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_bias_MSE(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                         DLT_scenario = DLT_scenario_full_description,
                                         EFF_scenario = EFF_scenario_full_description,
                                         omega = omega_num,
                                         agreement = agreement_text))),
          paste0("bias_MSE_parameters",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))

#optimal dose combination
write.csv(do.call("<-",list(paste0("optimal_dose_combination",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_optimal_dose_combination(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                         DLT_scenario = DLT_scenario_full_description,
                                                         EFF_scenario = EFF_scenario_full_description,
                                                         omega = omega_num,
                                                         agreement = agreement_text))),
          paste0("optimal_dose_combination",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))

#Power / type-I error

write.csv(do.call("<-",list(paste0("power_or_typeIerror",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_power(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                      type = power_or_typeIerror,
                                      DLT_scenario = DLT_scenario_full_description,
                                      EFF_scenario = EFF_scenario_full_description,
                                      omega = omega_num,
                                      agreement = agreement_text))),
          paste0("power_or_typeIerror",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))


#Stopping rules
write.csv(do.call("<-",list(paste0("stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_posterior_probabilities_stopping_rules(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                       delta_Z1 = delta_Z1,
                                                                       delta_Z2 = delta_Z2,
                                                                       delta_E0 = delta_E0,
                                                                       DLT_scenario = DLT_scenario_full_description,
                                                                       EFF_scenario = EFF_scenario_full_description,
                                                                       omega = omega_num,
                                                                       agreement = agreement_text))),
          paste0("stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#Sample size stopping rules
write.csv(do.call("<-",list(paste0("sample_size_stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_sample_size_stopping_rules(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                           DLT_scenario = DLT_scenario_full_description,
                                                           EFF_scenario = EFF_scenario_full_description,
                                                           omega = omega_num,
                                                           agreement = agreement_text))),
          paste0("sample_size_stopping_rules",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#Proportion patient optimal dose combinations with true P(pi_E) > theta_E
write.csv(do.call("<-",list(paste0("optimal_dose_combination_above_threshold",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_proportion_optimal_dose_combination_above_efficacy_threshold(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                                             DLT_scenario = DLT_scenario_full_description,
                                                                                             EFF_scenario = EFF_scenario_full_description,
                                                                                             omega = omega_num,
                                                                                             agreement = agreement_text))),
          paste0("optimal_dose_combination_above_threshold",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))




#Proportion patient allocation is dose combination with true P(pi_E) > theta_E
write.csv(do.call("<-",list(paste0("proportion_allocation",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_proportion_correct_patient_allocation(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                                      DLT_scenario = DLT_scenario_full_description,
                                                                      EFF_scenario = EFF_scenario_full_description,
                                                                      omega = omega_num,
                                                                      agreement = agreement_text))),
          paste0("proportion_allocation",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))


#Dose combinations to which patients were allocated
write.csv(do.call("<-",list(paste0("dose_combinations",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_dose_combinations(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                  DLT_scenario = DLT_scenario_full_description,
                                                  EFF_scenario = EFF_scenario_full_description,
                                                  omega = omega_num,
                                                  agreement = agreement_text))),
          paste0("dose_combinations",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))




# MTD curve
write.csv(do.call("<-",list(paste0("MTD_curve",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_MTD_curve(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                          DLT_scenario = DLT_scenario_full_description,
                                          EFF_scenario = EFF_scenario_full_description,
                                          omega = omega_num,
                                          agreement = agreement_text))),
          paste0("MTD_curve",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))



#DLT rate and proportion of trials with DLT rate above theta_Z + 0.1
write.csv(do.call("<-",list(paste0("DLT_rate",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis),
                            get_average_DLT_rate(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis))),
                                                 DLT_scenario = DLT_scenario_full_description,
                                                 EFF_scenario = EFF_scenario_full_description,
                                                 omega = omega_num,
                                                 agreement = agreement_text))),
          paste0("DLT_rate",DLT_scenario_variables,EFF_scenario_variables,agreement,omega,hypothesis,".csv"))

