## Functions necessary for the bayesian fitting algorithm

## Likelihood

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

# Default parameters
foc.param<-list(alphaC=0.05,alphaA=0.05,scaleConst=150,
                gamma0=0.0, gamma1=0.0,negReward0=0.0,negReward1=0.0,
                probFAA=1,interpReg=0,slopRegRelAC=0, slopRegPVL=0)


priors<-data.frame(best=as.numeric(foc.param),
                   lower=c(0,0,0,-1,-1,-50,-50,0,-50,-50,-50),
                   upper=c(1,1,500,1,1,50,50,1,50,50,50))

rownames(priors)<-names(foc.param)

## Priors
# Create own prior
densityPriorCreator <- function(par.sel){
  funcs<-list(
    alphaC=function(x) dbeta(x,2,5),
    alphaA=function(x) dbeta(x,2,5),
    scaleConst= function(x) dtruncnorm(x,0,500,150,100),
    gamma0= function(x) dtruncnorm(x,-1,1,0,0.2),# dbeta(pars[2],2,5)
    gamma1=function(x) dtruncnorm(x,-1,1,0,0.2),# dbeta(pars[2],2,5)
    negReward0=function(x) dtruncnorm(x,-50,50,0,15),
    negReward1=function(x) dtruncnorm(x,-50,50,0,15),
    probFAA=function(x) dbeta(x,2,2),
    interpReg= function(x)dtruncnorm(x,-50,50,0,15),
    slopRegRelAC=function(x)dtruncnorm(x,-50,50,0,15),
    slopRegPVL=function(x)dtruncnorm(x,-50,50,0,15)
  )
  funcs<-funcs[par.sel]
  priorFunc<-function(pars){
    do.call(sum,lapply(1:length(pars), function(x){
      funcs[[x]](pars[x])
    }))
  }
  return(priorFunc)
}

# Create uniform prior
# prior <- createUniformPrior(lower = priors[parSel,"lower"], 
#                             upper = priors[parSel,"upper"], 
#                             best = priors[parSel,"best"])


## Sampler for the prior distribution

samplePriorCreator <- function(par.sel){
  funcs<-list(
    alphaC=function(n) rbeta(n,2,5),
    alphaA=function(n) rbeta(n,2,5),
    scaleConst= function(n) rtruncnorm(n,0,500,150,100),
    gamma0= function(n) rtruncnorm(n,-1,1,0,0.2),# dbeta(pars[2],2,5)
    gamma1=function(n) rtruncnorm(n,-1,1,0,0.2),# dbeta(pars[2],2,5)
    negReward0=function(n) rtruncnorm(n,-50,50,0,15),
    negReward1=function(n) rtruncnorm(n,-50,50,0,15),
    probFAA=function(n) rbeta(n,2,5),
    interpReg= function(n) rtruncnorm(n,-50,50,0,15),
    slopRegRelAC=function(n) rtruncnorm(n,-50,50,0,15),
    slopRegPVL=function(n) rtruncnorm(n,-50,50,0,15)
  )
  funcs<-funcs[par.sel]
  sampleFunc<-function(n=1){
    sample.tmp<-do.call(cbind,lapply(1:length(funcs), function(x){
      funcs[[x]](n)
    }))
    colnames(sample.tmp)<-par.sel
    return(sample.tmp)
  }
  return(sampleFunc)
}

