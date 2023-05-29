library(here)
library(ggplot2)
source(here("loadData.R"))
# Cpp file with the simulation model
sourceCpp(here("ActCrit_R.cpp"))
source(here("bayes_funcs.R"))

logist2d<-function(x,y,pars){
  sapply(1:length(x), function(n)
    1/(1+exp(-(pars[1]+pars[2]*x[n]+pars[3]*y[n])))
  )
}

fieldData<-fread(here("Data","data_cleaner_abs_threa1.5.txt"))
names(fieldData)[4:8]<-c("abund_clean","abund_visitors","abund_resid",
                         "prob_Vis_Leav","group")


test.pars<-c(0,0,0)

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
  agentScen = 1,
  # Type of agent FAA (chuncking), PAA (not chuncking)
  propfullPrint = 0.7, 
  #Proportion of final rounds used to calculate predictions
  # alphaA,AlphaC, Gamma, NegRew, scalConst,probFAA
  # nRep=30, # number of replicate simulations
  Group = FALSE # grouped data according to social competence
  # from triki et al. 2020
)

test.data<-do_simulation(fieldData,foc.param,param_mcmc)

mean(test.data$agent)

logist2d(test.data$rel_abund_clean,test.data$prob_Vis_Leav,test.pars)
