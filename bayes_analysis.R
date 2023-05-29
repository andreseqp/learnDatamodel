## Analize and visualize MCMC chains

library(here)
library(BayesianTools)
library(truncnorm)
library(Rcpp)
library(data.table)
sourceCpp(here("ActCrit_R.cpp"))
source(here("bayes_funcs.R"))

defaultPars<-foc.param
scenarios_select<-list(BT_gam_Nrew_sca=c("scaleConst", "gamma0","negReward0"),
                       BT_gam_sca=c("scaleConst", "gamma0"),
                       BT_Nrew_sca=c("scaleConst", "negReward0"))
fieldData<-fread(here("Data","data_cleaner_abs_threa1.5.txt"))
names(fieldData)[4:8]<-c("abund_clean","abund_visitors","abund_resid",
                         "prob_Vis_Leav","group")
param_mcmc<-list(
  totRounds=10000, 
  # Number of rounds in the learning model
  ResReward=1,VisReward=1, 
  # Magnitude of reward for residents and visitors
  ResProbLeav=0, 
  # Prob. of resident leaving the station if unattended
  scenario=0, 
  # scenario of how the clients reach the station
  #nature 0, experiment 1, marketExperiment 2, ExtendedMarket 3
  inbr=0,outbr=0, 
  # probability of clients seeking similar/ different,
  # clients, respectively, in the station
  seed=1,  
  # Seed for the random number generator
  agentScen = 0,
  # Type of agent FAA (chuncking), PAA (not chuncking)
  propfullPrint = 0.7, 
  #Proportion of final rounds used to calculate predictions
  # alphaA,AlphaC, Gamma, NegRew, scalConst,probFAA
  nRep=30, # number of replicate simulations
  Group = FALSE, # grouped data according to social competence
  # from triki et al. 2020
)

# FAA --------------------------------------------------------------------------

scenario <-"BT_Nrew_sca"
parSel <- scenarios_select[[scenario]]
# Load MCMC chains
MCMC.FAA.loaded <- readRDS(file = here(paste0(scenario,"_"),
                                     "MCMC_FAA.rda"))

# Plot the MCMC chains ----------------------------------------------------------

par()
plot(MCMC.FAA.loaded)
summary(MCMC.FAA.loaded)
marginalPlot(MCMC.FAA.loaded)
gelmanDiagnostics(MCMC.FAA.loaded)

# PAA --------------------------------------------------------------------------

scenario <-"BT_PAA_gam_Nrew_sca"
# Load MCMC chains
MCMC.PAA.loaded <- readRDS(file = here(paste0(scenario,"_"),
                                       "MCMC_FAA.rda"))

# Plot the MCMC chains ----------------------------------------------------------

par()
plot(MCMC.PAA.loaded)
summary(MCMC.PAA.loaded)
marginalPlot(MCMC.PAA.loaded)
gelmanDiagnostics(MCMC.PAA.loaded)
