# -------------- generate json files with parameters------------------------------------------ #

library("jsonlite")
library("here")
source(here("loadData.R"))

simsDir<-"Simulations"

clusterHome<-"/home/ubuntu/learnData"

fileName<-"parameters.json"

# name of the scenario for which parameter files will be produces
# format used for the name: MCMCclean_*parameter1_paramerter2_parameter3.....
scenario<- "MCMCclean_PAA2_gam_Nrew_sca"#""


# MCMC fit - Generate json parameter files for -------------------------------------

# for MCMC
param_mcmc<-list(totRounds=10000, 
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
                 # Type of agent FAA (chuncking), PAA (not chuncking)
                 propfullPrint = 0.7, 
                 #Proportion of final rounds used to calculate predictions
                 sdPert=c(0.05,0.05,0.2,5,200), 
                 # Width of the perturbation kernel for each parameter
                 # alphaA,AlphaC, Gamma, NegRew, scalConst
                 chain_length=100000, # Chain length
                 init=c(0.05,0.05,0,0,30), # Initial values for each of the parameters
                 # alphaA,AlphaC, Gamma, NegRew, scalConst
                 pertScen = c(FALSE,FALSE,TRUE,TRUE,TRUE), # boolean controlling which
                 # parameter is perturbed
                 MCMC =1, # run chain=1, run prediction for one parameter set = 0
                 nRep=1, # number of replicate simulations
                 Group = FALSE, # grouped data according to social competence
                 # from triki et al. 2020
                 dataFile =
                 paste(clusterHome,"Data","data_cleaner_abs_threa1.5.txt",
                                  sep="/"),
                   # here("Data","data_cleaner_abs_threa1.5.txt"),
                 # location of the data file
                 # here("Simulations",
                 #                 paste0(scenario,"_"),"data_MCMC_fake.txt"))
                 folderL=paste0(here(simsDir),"/",scenario,"_/")
                 # Local folder
                 )

# Make folder
check_create.dir(here(simsDir),param = rep(scenario,1),
                 values = c(""))


# Loop to make parameter files ------------------------------------------------
# and give different starting values

for(seed in 5:5){
  param_mcmc$folder<-param_mcmc$folderL
  param_mcmc$folder<-paste(clusterHome,paste0(scenario,"_/"),sep="/")
  param_mcmc$init<-c(0.05,0.05,runif(1,max = 0.6),runif(1,max = 1,min = 0),
                     runif(1,max = 100,min = 5))
  param_mcmc$init[3:4]<-param_mcmc$init[3:4]*param_mcmc$pertScen[3:4]
  param_mcmc$seed <- seed
  fileName<-paste("parametersMCMCcluster_",seed,".json",sep="")
  outParam<-toJSON(param_mcmc,auto_unbox = TRUE,pretty = TRUE)
  if(file.exists(paste(param_mcmc$folderL,fileName,sep = ''))){
    currFile<-fromJSON(paste(param_mcmc$folderL,fileName,sep = ''))
    if(sum(unlist(currFile)!=unlist(param_mcmc))>0){
      write(outParam,paste(param_mcmc$folderL,fileName,sep = "/"))
      
    }
  }else{
    write(outParam,paste(param_mcmc$folderL,fileName,sep = ""))
  }
}

# Predicted values -------------------------------------------------------------

# Run the learning model with parameters estimated through the MCMC

## Using the mode of the marginal posterior distributions

scenario<-"MCMCclean_nRew_sca"


# Load MCMC data
MCMCdata<-loadMCMCrun(paste0(scenario,"_"))
# get the posterior of gamma
densGamma<-lapply(c("gamma"),function(x)density(MCMCdata[,get(x)]))
names(densGamma)<-c("gamma")
# compute the mode of the posterior
modeGamma<-lapply(c(1),function(id)
  densGamma[[id]]$x[densGamma[[id]]$y==max(densGamma[[id]]$y)])
names(modeGamma)<-c("gamma")
# get the posterior of eta
densNR<-lapply(c("negReward"),function(x)density(MCMCdata[,get(x)]))
names(densNR)<-c("NR")
modeNegrew<-lapply(c(1),function(id)
  densNR[[id]]$x[densNR[[id]]$y==max(densNR[[id]]$y)])
names(modeNegrew)<-c("NR")
# get the posterior of scaling constant
densScal<-density(MCMCdata$scaleConst)
modescal<-densScal$x[densScal$y==max(densScal$y)]

# parameter file 
param_pred<-list(totRounds=10000,ResReward=1,VisReward=1,
            ResProbLeav=0,scenario=0, inbr=0,outbr=0,forRat=0.0,
            seed=1, propfullPrint = 0.7,sdPert=c(0.05,0.05,0.1,0.05,1),
            chain_length=0,
            init=c(0.05,0.05,0,0,35),# alphaA,AlphaC, Gamma, NegRew, scalConst
            init2=c(0.05,0.05,0,0,35),
            pertScen = c(FALSE,FALSE,FALSE,TRUE,TRUE), 	
            MCMC =0, 
            nRep=30, # number of replicates to run
            dataFile = here("Data","data_cleaner_abs_threa1.5.txt"),
            Group = FALSE,
            folderL=paste0(here(simsDir),"/",scenario,"_/"))



# create dir
check_create.dir(here(simsDir,paste0(scenario,"_")),param = rep("samplesPost",1),
                 values = c(""))

# Change is running on cluster
param_pred$folder<-param_pred$folderL

# For samples from the posterior
param_pred_samp<-param_pred
param_pred_samp$folderL<-paste0(here(simsDir),"/",scenario,"_/","samplesPost_/")
param_pred_samp$folder<-param_pred_samp$folderL

## Use modes for predictions
param_pred$init[c(3,4,5)]<- c(modeGamma$gamma,modeNegrew$NR,modescal)
                            

nsamples<-100

# Get samples
postSamp<-MCMCdata[sample(x = 1:dim(MCMCdata)[1],
                          size = nsamples,replace = FALSE),]

# Print file using the mode
fileName<-"parameters_pred_1.json"
outParam.pred<-toJSON(param_prep,auto_unbox = TRUE,pretty = TRUE)
write(outParam.pred,paste(param_pred$folderL,fileName,sep = ""))


# print file using the samples
for(i in 1:nsamples){
  param_pred_samp$init[c(3,4,5)]<- as.numeric(postSamp[i,
                                        .(gamma,negReward,scaleConst)])
  param_pred_samp$init[3:4]<-param_pred_samp$init[3:4]*param_pred_samp$pertScen[3:4]
  # param_pred$seed<- i
  fileName<-paste("parameters_pred_",i,".json",sep="")
  outParam.pred<-toJSON(param_pred_samp,auto_unbox = TRUE,pretty = TRUE)
  if(file.exists(paste(param_pred_samp$folderL,fileName,sep = ''))){
    currFile<-fromJSON(paste(param_pred_samp$folderL,fileName,sep = ''))
    if(sum(unlist(currFile)!=unlist(param_pred_samp))>0){
      write(outParam.pred,paste(param_pred_samp$folderL,fileName,sep = "/"))
    }
  }else{
    write(outParam.pred,paste(param_pred_samp$folderL,fileName,sep = ""))
  }
}

# For the countour plots ------------------------------------------------------
# To generate the contour plots, predictions need to b obtained varying 
# systematically the values of cleaner abundance and visitor leaving probability

param<-list(totRounds=10000,ResReward=1,VisReward=1,
            ResProb=c(0.2),
            VisProb=c(0.2),
            ResProbLeav=0,VisProbLeav=1,negativeRew=-modeNegrew$NR,
            scenario=0,
            inbr=0,outbr=0,trainingRep=10,forRat=0.0,
            alphaT=0.05,printGen=1,seed=1, gammaRange=I(c(modeGamma$gamma)),
            netaRange=I(c(1)),alphaThRange=I(c(0.05)),numlearn=1,
            propfullPrint = 0.7,
            alphaThNch=0.05,
            folderL=paste0(here(simsDir),scenario,"_/"))

# Arrays with the values of cleaner abundance and visitor leaving probability
rangLeav<-seq(0.02,0.4,length.out = 10)
rangAbund<-seq(0.05,0.9,length=10)
# rangScen<-c(0)
# rangAlphNC<-c(0,0.5,1)

# General folder for analysis
check_create.dir(here(simsDir),param = rep(scenario,1),
                 values = c(""))

listfolders<-check_create.dir(here("Simulations",paste0(scenario,"_")),
                              param = rep("Vlp",length(rangLeav)),
                              values = round(rangLeav,2))



# Loop through parameter names and values creating JSONs -----------------------

for (i in 1:length(rangLeav)) {
  for(j in 1:length(rangAbund)){
    param$negativeRew<-param$negativeRew*param_pred$pertScen[4]
    param$gammaRange<-I(c(param$gammaRange*param_pred$pertScen[3]))
    param$folderL<-paste0(here("Simulations",paste0(scenario,"_"),listfolders[i]),"/")
    param$folder<-param$folderL 
    param$ResProb<-c((1-rangAbund[j])/2)
    param$VisProb<-c((1-rangAbund[j])/2)
    param$VisProbLeav<-rangLeav[i]
    outParam<-toJSON(param,auto_unbox = TRUE,pretty = TRUE)
    fileName<-paste("parameters_contour_",j,".json",sep="")
    if(file.exists(paste(param$folderL,fileName,sep = ''))){
      currFile<-fromJSON(paste(param$folderL,fileName,sep = ''))
      if(sum(unlist(currFile)!=unlist(param))>0){
        write(outParam,paste(param$folderL,fileName,sep = "/"))
      }
    }
    else{
      write(outParam,paste(param$folderL,fileName,sep = ""))
    }
    
  }
}# 

