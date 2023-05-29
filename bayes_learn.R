# install.packages(c(
#   "Rcpp","here","BayesianTools","data.table","jsonlite","truncnorm"
# ))
# suppressPackageStartupMessages(c(
#   library(Rcpp,lib.loc="RLibs"),
#   library(here,lib.loc="RLibs"),
#   library("BayesianTools",lib.loc="RLibs"),
#   library(data.table,lib.loc="RLibs"),
#   library("jsonlite",lib.loc="RLibs"),
#   # library("RcppJson")
#   library("ps",lib.loc="RLibs"),
#   library(devtools,lib.loc="RLibs"),
#   library("truncnorm",lib.loc="RLibs")
# ))

suppressPackageStartupMessages(c(
  library(Rcpp),
  library(here),
  library("BayesianTools"),
  library(data.table),
  library("jsonlite"),
  # library("RcppJson")
  library("ps"),
  library(devtools),
  library("truncnorm")
))

source(here("loadData.R"))
# Cpp file with the simulation model
sourceCpp(here("ActCrit_R.cpp"))

nIter<-100
nChainsI<-1

#load data
fieldData<-fread(here("Data","data_cleaner_abs_threa1.5.txt"))
names(fieldData)[4:8]<-c("abund_clean","abund_visitors","abund_resid",
                         "prob_Vis_Leav","group")


# choosing which parameters to calibrate
# parSel = c("scaleConst", "gamma0","negReward0")

scenarios_select<-list(BT_log_sca=
                         c("scaleConst", "interpReg","slopRegRelAC",
                           "slopRegPVL"))
# ,
#                 BT_log_Nrew_sca=c("scaleConst", "negReward0"))



scenarios<-names(scenarios_select)

# for(scenario in scenarios){
  
  scenario<-scenarios[1]
  parSel <- scenarios_select[[scenario]]

  check_create.dir(here(),param = rep(scenario,1),
                 values = c(""))
  
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
  
  
  source(here("bayes_funcs.R"))
  
  defaultPars<-foc.param
  
  # define priors and their sampling function
  densityPrior<-densityPriorCreator(parSel)
  samplerPrior<-samplePriorCreator(parSel)
  prior <- createPrior(density = densityPrior, sampler = samplerPrior,
                       lower = priors[parSel,"lower"], 
                       upper = priors[parSel,"upper"], 
                       best = NULL)
  
  
  # library(ggplot2)
  # library(dplyr)
  # 
  # prior.short<-samplerPrior(10000) %>% as.data.table()
  # 
  # prior.short %>% melt() %>%
  #   ggplot(aes(x=value,fill=variable)) +
  #   geom_histogram() +
  #   facet_grid(~variable,scales = "free_x")+
  #   theme_classic()
  
  
  
## Plotting the priors ---------------------------------------------------------
# par(mfrow=c(3,1))
# plot(x=seq(0,500,by=0.1),
#      y=dtruncnorm(x=seq(0,500,by=0.1),0,500,150,100),type="l")
# 
# plot(x=seq(-1,1,by=0.01),
#      y=dtruncnorm(seq(-1,1,by=0.01),-1,1,0,0.2),type="l")
# 
# plot(x=seq(-50,50,by=0.01),
#      y=dtruncnorm(seq(-50,50,by=0.01),-50,50,0,15),type="l")


# Set up the bayesian engine ---------------------------------------------------
  bayesianSetup <- createBayesianSetup(LogLihood, prior, 
                                     names = parSel)

# settings for the sampler, iterations should be increased for real applicatoin
  settings <- list(iterations = nIter, nrChains = 1)

  ## Write out file with parameter values 
  fileName<-paste("parametersMCMC_",".json",sep="")
  outParam<-toJSON(param_mcmc,auto_unbox = TRUE,pretty = TRUE)
  write(outParam,here(paste0(scenario,"_"),fileName))
  

  # sourceCpp(here("ActCrit_R.cpp"))
  do_simulation(fieldData,foc.param,param_mcmc)
  settings <- list(iterations = 100, nrChains = 1)
  tmp.mcmc<-runMCMC(bayesianSetup, settings, sampler = "DEzs")

  

  ## Run MCMC chains in parallel
  ## Start cluster with n cores for n chains and export BayesianTools library
  nChains<-nChainsI
  cl <- parallel::makeCluster(nChains+1,outfile="out")
  parallel::clusterExport(cl,c("param_mcmc",
                               "foc.param",
                               "parSel",
                               "fieldData")
                          )
  parallel::clusterEvalQ(cl, {
    # library(BayesianTools,lib.loc="RLibs")
    # library(Rcpp,lib.loc="RLibs")
    # library(here,lib.loc="RLibs")
    # library("truncnorm",lib.loc="RLibs")
    # library(data.table,lib.loc="RLibs")
    # library("jsonlite",lib.loc="RLibs")
    
    library(BayesianTools)
    library(Rcpp)
    library(here)
    library("truncnorm")
    library(data.table)
    library("jsonlite")
    
    # library("RcppJson",lib.loc="RLibs")
    # Cpp file with the simulation model
    sourceCpp(here("ActCrit_R.cpp"))
    
    defaultPars<-foc.param
    }
  )

  ## calculate parallel n chains, for each chain the likelihood will be calculated on one core
  MCMC.FAA <- parallel::parLapply(cl, 1:nChains, 
                                  fun = function(X, bayesianSetup, settings) 
    runMCMC(bayesianSetup, settings, sampler = "DEzs"), bayesianSetup, settings)
  
  
  ## Combine the chains
  MCMC.FAA <- createMcmcSamplerList(MCMC.FAA)

  # Save files for future analysis
  saveRDS(MCMC.FAA, file= here(paste0(scenario,"_"),"MCMC_FAA.rda"))

  print(paste0("Done ",scenario))
}

