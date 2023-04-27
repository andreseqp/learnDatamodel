library(Rcpp)
library(here)
library("BayesianTools")
library(data.table)
library("jsonlite")
library("RcppJson")
library("truncnorm")
# Cpp file with the simulation model
sourceCpp(here("ActCrit_R.cpp"))


fieldData<-fread(here("Data","data_cleaner_abs_threa1.5.txt"))
names(fieldData)[4:8]<-c("abund_clean","abund_visitors","abund_resid",
                         "prob_Vis_Leav","group")

# Parameters not to be fitted
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
   forRat=0.0, 
   # Rate at which cleaners forget what they have learned
   seed=1,  
   # Seed for the random number generator
   agent="FAA",
   agentScen = 0,
   # Type of agent FAA (chuncking), PAA (not chuncking)
   propfullPrint = 0.7, 
   #Proportion of final rounds used to calculate predictions
   # alphaA,AlphaC, Gamma, NegRew, scalConst,probFAA
   nRep=30, # number of replicate simulations
   Group = FALSE, # grouped data according to social competence
   # from triki et al. 2020
   groupPars = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
)

foc.param<-list(alphaC=0.05,alphaA=0.05,scaleConst=150,
                gamma0=0.9, gamma1=0.9,negReward0=0.0,negReward1=0.0,
                probFAA=1,interpReg=1,slopRegRelAC=2,
                slopRegPVL=2)

priors<-data.frame(best=as.numeric(foc.param),
                   lower=c(0,0,0,0,0,-50,-50,0,-20,-20,-20),
                   upper=c(1,1,500,1,1,50,50,1,20,20,20))

rownames(priors)<-names(foc.param)

defaultPars<-foc.param

# choosing which parameters to calibrate
parSel = c("scaleConst", "gamma0", "negReward0")


LogLihood<-function(pars){
    # set parameters that are not calibrated on default values 
    x = defaultPars
    x[parSel] = pars[parSel]
    predicted <- do_simulation(emp_data = fieldData,
                  focal_param = x,
                  sim_param =   param_mcmc)    # run simulations
    logliHood<-sum(dbinom(predicted[,"score_visitor"],20,
           predicted[,"marketPred"],log = TRUE)) # 
    # Log likelihood calculated using the number of visitor choices, 
    # the total number of trials (20) and the model prediction
    return(logliHood)
}

# Create uniform prior
prior <- createUniformPrior(lower = priors[parSel,"lower"], 
                            upper = priors[parSel,"upper"], 
                            best = priors[parSel,"best"])

# Create own prior
densityPrior <- function(pars){
  scaleConstD <-dtruncnorm(pars[1],0,500,150,100)
  gamma1D <- dbeta(pars[2],2,5)
  negreward1D <-dtruncnorm(pars[3],-50,50,0,15)
  return(scaleConstD+gamma1D+negreward1D)
}

samplerPrior <- function(n=1){
  scaleConstR <-rtruncnorm(n,0,500,150,100)
  gamma1R <- rbeta(n,2,5)
  negreward1R <-rtruncnorm(n,-50,50,0,15)
  return(cbind(scaleConstR,gamma1R,negreward1R))
}

prior <- createPrior(density = densityPrior, sampler = samplerPrior,
                     lower = priors[parSel,"lower"], 
                     upper = priors[parSel,"upper"], 
                     best = NULL)



# Set up the bayesian engine
bayesianSetup <- createBayesianSetup(LogLihood, prior, 
                                     names = parSel)

# settings for the sampler, iterations should be increased for real applicatoin
settings <- list(iterations = 100000, nrChains = 1)

MCMC.FAA <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)


par()
plot(MCMC.FAA)
summary(MCMC.FAA)
marginalPlot(MCMC.FAA)
gelmanDiagnostics(MCMC.FAA)

## Run MCMC chains in parallel

## Start cluster with n cores for n chains and export BayesianTools library
nChains<-5
cl <- parallel::makeCluster(nChains)
parallel::clusterEvalQ(cl, {
  library(BayesianTools)
  library(Rcpp)
  library(here)
  library("truncnorm")
  library(data.table)
  library("jsonlite")
  library("RcppJson")
  # Cpp file with the simulation model
  sourceCpp(here("ActCrit_R.cpp"))
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
    forRat=0.0, 
    # Rate at which cleaners forget what they have learned
    seed=1,  
    # Seed for the random number generator
    agent="FAA",
    agentScen = 0,
    # Type of agent FAA (chuncking), PAA (not chuncking)
    propfullPrint = 0.7, 
    #Proportion of final rounds used to calculate predictions
    # alphaA,AlphaC, Gamma, NegRew, scalConst,probFAA
    nRep=30, # number of replicate simulations
    Group = FALSE, # grouped data according to social competence
    # from triki et al. 2020
    groupPars = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
  )
  defaultPars<-list(alphaC=0.01,alphaA=0.01,scaleConst=150,
                    gamma0=0.9, gamma1=0.9,negReward0=0.0,negReward1=0.0,
                    probFAA=1,interpReg=1,slopRegRelAC=2,
                    slopRegPVL=2)
  parSel = c("scaleConst", "gamma0", "negReward0")
  fieldData<-fread(here("Data","data_cleaner_abs_threa1.5.txt"))
  names(fieldData)[4:8]<-c("abund_clean","abund_visitors","abund_resid",
                           "prob_Vis_Leav","group")
}
)

## calculate parallel n chains, for each chain the likelihood will be calculated on one core
MCMC.FAA <- parallel::parLapply(cl, 1:nChains, fun = function(X, bayesianSetup, settings) 
  runMCMC(bayesianSetup, settings, sampler = "DEzs"), bayesianSetup, settings)

## Combine the chains
MCMC.FAA <- createMcmcSamplerList(MCMC.FAA)

head(MCMC.FAA)


par()
plot(MCMC.FAA)
summary(MCMC.FAA)
marginalPlot(MCMC.FAA)
gelmanDiagnostics(MCMC.FAA)

