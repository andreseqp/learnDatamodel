library(Rcpp)
library(here)
library("BayesianTools")
library(data.table)
library("jsonlite")
library("RcppJson")


foc.param<-list(alphaC=0.01,alphaA=0.01,scaleConst=150,
     gamma=c(0.1,0.2),negReward=c(0.5,0.8),
     probFAA=c(1,1),interpReg=1,slopRegRelAC=2,
     slopRegPVL=2)

fieldData<-fread(here("Data","data_cleaner_abs_threa1.5.txt"))
names(fieldData)[4:8]<-c("abund_clean","abund_visitors","abund_resid",
                         "prob_Vis_Leav","group")

scenario<-"testBayesianTools"

param_mcmc<-list(
    totRounds=2, 
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
   agent="PAA",
   agentScen = 0,
   # Type of agent FAA (chuncking), PAA (not chuncking)
   propfullPrint = 0.7, 
   #Proportion of final rounds used to calculate predictions
   sdPert=c(0.05,0.05,0.25,6.5,100,0.2), 
   # Width of the perturbation kernel for each parameter
   # alphaA,AlphaC, Gamma, NegRew, scalConst,probFAA
   chain_length=100000, # Chain length
   init=c(0.05,0.05,0,0,30,0), # Initial values for each of the parameters
   # alphaA,AlphaC, Gamma, NegRew, scalConst,probFAA
   pertScen = c(FALSE,FALSE,FALSE,TRUE,TRUE,FALSE), # boolean controlling which
   # parameter is perturbed
   MCMC =1, # run chain=1, run prediction for one parameter set = 0
   nRep=1, # number of replicate simulations
   Group = FALSE, # grouped data according to social competence
   # from triki et al. 2020
   groupPars = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
   dataFile =
   # paste(clusterHome,"Data","data_cleaner_abs_threa1.5.txt",
   #                  sep="/"),
   here("Data","data_cleaner_abs_threa1.5.txt"),
   # location of the data file
   # here("Simulations",
   #                 paste0(scenario,"_"),"data_MCMC_fake.txt"))
   folderL=paste0(here("Simulations"),"/",scenario,"_/")
   # Local folder
)


param_mcmc$agentScen <-0


sourceCpp(here("SWE_ActCrit.cpp"))

do_simulation(emp_data = fieldData,
              focal_param = foc.param,
              # fileStr= outParam)
              sim_param =   param_mcmc)


sourceCpp(here("test.cpp"))

do_simulation(emp_data = fieldData,
              focal_param = foc.param,
              # fileStr= outParam)
              sim_param =   param_mcmc)

sourceCpp(here("ActCrit_R.cpp"))

do_simulation_test(emp_data = fieldData,
                   focal_param = foc.param,
                   sim_param =   param_mcmc)


predicton<-do_simulation(emp_data = fieldData,
              focal_param = foc.param,
              # fileStr= outParam)
              sim_param =   param_mcmc)

timesTwo(3)

?VSEM

focal<-new(model_param)
focal$set_gamma(0.5,0.5)
focal$logist()
focal$methodJson()

showPi(focal$myJson)

new_param<-new(model_param)
new_param$set_gamma(1,10)
new_param$logist()

focal$copy(new_param$get_ptr())

focal$logist()


focal$set_gamma(1,1)
focal$logist()
focal$set_gamma(0.5,0.6)
focal$logist()

hello()
bla()
bla2(42, 0.42)

w <- new(World)
w$greet()
w$set("hohoho")
w$greet()

muUnif<-new(Uniform,0,2)

muUnif$draw(5)
