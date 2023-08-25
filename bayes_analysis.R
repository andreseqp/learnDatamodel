## Analize and visualize MCMC chains

library(here)
library(BayesianTools)
library(truncnorm)
library(Rcpp)
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
sourceCpp(here("ActCrit_R.cpp"))
source(here("bayes_funcs.R"))

defaultPars<-foc.param
scenarios_select<-list(BT_FAA_gam_Nrew_sca=
                         c("scaleConst", "gamma0","negReward0"),
                       BT_FAA_gam_sca=c("scaleConst", "gamma0"),
                       BT_FAA_Nrew_sca=c("scaleConst", "negReward0"))
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
  Group = FALSE # grouped data according to social competence
  # from triki et al. 2020
)

# FAA --------------------------------------------------------------------------

scenario <-"BT_FAA_gam_sca"
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

scenario <-"BT_FAA_gam_sca"
# Load MCMC chains
MCMC.PAA.loaded <- readRDS(file = here(paste0(scenario,"_"),
                                       "MCMC_FAA.rda"))

# Plot the MCMC chains ----------------------------------------------------------

par()
plot(MCMC.PAA.loaded)
summary(MCMC.PAA.loaded)
marginalPlot(MCMC.PAA.loaded)
gelmanDiagnostics(MCMC.PAA.loaded)

# Compare two models via Bayes factor ------------------------------------------

ML.FAA <- marginalLikelihood(MCMC.FAA.loaded)
ML.PAA <- marginalLikelihood(MCMC.PAA.loaded)

exp(M1$ln.ML) / ( exp(M1$ln.ML) + exp(M2$ln.ML))

# Why is there not likelihood signal?

test.param <- foc.param
test.param$scaleConst<-250
test.param$alphaC<-0.05
test.param$alphaA<-0.05
param_mcmc.1<-param_mcmc
param_mcmc.1$totRounds<-20000

mean(fieldData$score_visitor/20)

sum(dbinom(fieldData$score_visitor,20,
           fieldData$score_visitor/20,log = TRUE))

sourceCpp(here("ActCrit_R.cpp"))


parSel<- names(test.param)

hist.data<-ggplot(fieldData,aes(x=score_visitor/20))+geom_histogram()+
  xlim(0,1)+ggtitle(paste0("Log-likelihood_data: ",
                           sum(dbinom(fieldData$score_visitor,20,
                             fieldData$score_visitor/20,log = TRUE))))


test.param.0<-test.param
test.param.0$gamma0<-0.3
test.param.0$probFAA<-0

FAA_gam_0 <- do_simulation(fieldData,test.param.0,param_mcmc.1)
mean(FAA_gam_0[,"marketPred"])
hist.0<-ggplot(FAA_gam_0,aes(x=marketPred))+geom_histogram()+
  ggtitle(paste0("Log-likelihood_0: ",LogLihood(pars = test.param.0)))+
  xlim(0,1)

test.param.1<-test.param
test.param.1$gamma0<-0.5
test.param.1$probFAA<-1
FAA_gam_1 <- do_simulation(fieldData,test.param.1,param_mcmc.1)
mean(FAA_gam_1[,"marketPred"])
LogLihood(pars = test.param.1)
hist.1<-ggplot(FAA_gam_1,aes(x=marketPred))+geom_histogram()+
  ggtitle(paste0("Log-likelihood_1: ",LogLihood(pars = test.param.1)))+
  xlim(0,1)
test.param.2<-test.param
test.param.2$gamma0<-0
test.param.2$negReward0<-0.5
test.param.2$probFAA<-1
FAA_gam_2 <- do_simulation(fieldData,test.param.2,param_mcmc.1)
mean(FAA_gam_2[,"marketPred"])
LogLihood(pars = test.param.2)
hist.2<-ggplot(FAA_gam_2,aes(x=marketPred))+geom_histogram()+
  ggtitle(paste0("Log-likelihood_2: ",LogLihood(pars = test.param.2)))+
  xlim(0,1)
test.param.3<-test.param
test.param.3$gamma0<-0
test.param.3$negReward0<-0.5
test.param.3$probFAA<-0
FAA_gam_3 <- do_simulation(fieldData,test.param.3,param_mcmc.1)
mean(FAA_gam_3[,"marketPred"])
LogLihood(pars = test.param.3)
hist.3<-ggplot(FAA_gam_3,aes(x=marketPred))+geom_histogram()+
  ggtitle(paste0("Log-likelihood_3: ",LogLihood(pars = test.param.3)))+
  xlim(0,1)
test.param.4<-test.param
test.param.4$gamma0<-0.5
test.param.4$negReward0<-0.5
test.param.4$probFAA<-1
FAA_gam_4 <- do_simulation(fieldData,test.param.4,param_mcmc.1)
mean(FAA_gam_4[,"marketPred"])
LogLihood(pars = test.param.4)
hist.4<-ggplot(FAA_gam_4,aes(x=marketPred))+geom_histogram()+
  ggtitle(paste0("Log-likelihood_4: ",LogLihood(pars = test.param.4)))+
  xlim(0,1)

plot_grid(hist.data,hist.0,hist.1,hist.2,
          hist.3,hist.4)




str(FAA_gam_0)
test.param$probFAA<-1
test.param$gamma0<-1
test.param$gamma1<-0.9
FAA_gam_1 <- do_simulation(fieldData,test.param,param_mcmc)
mean(FAA_gam_1[,"marketPred"])

test.param$probFAA<-1
test.param$gamma0<-1
FAA_gam_2 <- do_simulation(fieldData,test.param,param_mcmc)
mean(FAA_gam_2[,"marketPred"])

param_mcmc.1<-param_mcmc
param_mcmc.1$agentScen<-1
test.param$gamma0<-1
FAA_gam_2 <- do_simulation(fieldData,test.param,param_mcmc.1)
mean(FAA_gam_2[,"marketPred"])



param_mcmc$agentScen<-1
test.param$gamma0<-0.9
FAA_gam_2 <- do_simulation(fieldData,test.param,param_mcmc)



plot(x=FAA_gam_0$marketPred,y=FAA_gam_2$marketPred)

par(mfrow=x(2,1))
plot(marketPred~scodata=FAA_gam_0)

FAA_gam_0<-do_simulation(fieldData,test.param,param_mcmc)
