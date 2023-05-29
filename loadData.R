# ------------------------ Load data ------------------------------------------#

# Libraries

require("data.table")
require("jsonlite")
require("rlist")

# get files --------------------------------------------------------------------

getFilelist<-
  # reads de list of files and filters it according to a list of parameters
  # and values of interest
  function(folder, # folder where the files are
           listparam=NULL, 
           # list strings providing the parameters 
           # of interest
           values=NULL, 
           # list of values matching the list in 
           # listparam
           fullNam=FALSE # Boolean. Append full path name
           ){
  posAgen<-c("PAA","FAA","DP","p1")
  listRaw<-list.files(folder,recursive = TRUE,full.names = fullNam)
  if(length(listRaw)<1)warning("No files in such location",immediate. = TRUE)
  listRaw<-grep(".txt",listRaw,value = TRUE)
  fullList<-vector("list",4)
  names(fullList)<-posAgen
  if(length(listparam)!=length(values)){
    # parameter and value lists must be the same length
    warning("Parameter list and values don't match",immediate. = TRUE)
  }
  else{
      if(is.null(listparam)){
        # if there is no list, use all the data
        paramList<-listRaw
  }
    else{
      # if there is list, filter the data
      #regExpList<-paste0(listparam,values,"_",sep="")
      for (param in unique(listparam)){
        valsparam<-values[grep(param,listparam)]
        listRaw<-do.call(list.append,lapply(paste(param,valsparam,"_",sep=""),
                                           grep,x=listRaw,value=TRUE))
      }
      
      paramList<-listRaw
    }
      for(agent in posAgen){
        # Create a list with the lists separated by type of agents
        listAgent<-grep(agent,paramList,value = TRUE)
        fullList[[agent]]<-listAgent
      }
      return(fullList)
    }
}


# Load raw data ----------------------------------------------------------------

loadRawData<-function(folder,agent,listparam,values){
  setwd(folder)
  fullList<-getFilelist(folder,listparam,values)
  DT<-do.call(rbind,lapply(fullList[[agent]],fread))
  DT$option<-ifelse((DT$Client1==1 & DT$Client2==0) | 
                   (DT$Client1==0 & DT$Client2==1),"RV",NA)
  DT$option<-ifelse((DT$Client1==0 & DT$Client2==0),"RR",DT$option)
  DT$option<-ifelse((DT$Client1==1 & DT$Client2==1),"VV",DT$option)
  DT$option<-ifelse((DT$Client1==0 & DT$Client2==2) | 
                   (DT$Client1==2 & DT$Client2==0),"R0",DT$option)
  DT$option<-ifelse((DT$Client1==1 & DT$Client2==2) | 
                   (DT$Client1==2 & DT$Client2==1),"V0",DT$option)
  DT$option<-ifelse((DT$Client1==2 & DT$Client2==2),"00",DT$option)
  
  return(DT)
}

# Load the JSON files with parameter values  -----------------------------------

getParam<-function(folder,agent,listparam=NULL,values=NULL){
  setwd(folder)
  irrelPar<-c("gamma","tau","neta","pR","pV")
  listRaw<-list.files(folder,recursive = TRUE)
  jsonsList<-grep(".json",listRaw,value = TRUE)
  indRelPar<-seq(length(listparam))
  for(param in irrelPar){
      indRelPar<-grep(param,listparam,invert = TRUE)
      listparam<-listparam[indRelPar]
      values<-values[indRelPar]
  }
  if(length(listparam)!=length(values)){
    warning("Parameter list and values don't match",immediate. = TRUE)
  }
  else{
    if(is.null(listparam)){
      finalList<-jsonsList
    }
    else{
      for (param in unique(listparam)){
        valsparam<-values[grep(param,listparam)]
        jsonsList<-do.call(list.append,lapply(paste(param,valsparam,"_",sep=""),
                                            grep,x=jsonsList,value=TRUE))
      }
    }
  }
  jsons<-do.call(list,lapply(jsonsList,fromJSON))
  return(jsons)
}


# Load time intervals ----------------------------------------------------------
# divides raw data in interval according to interV and calculates 
# the RV preference is such interval

file2timeInter<-function(filename,interV,maxAge=-2){
  extPar<-strsplit(filename,split ="_/")[[1]][1]
  parVal<-as.numeric(gsub("[[:alpha:]]",extPar,replacement = ''))
  extPar<-gsub("[[:digit:]]",extPar,replacement = '')
  tmp<-fread(filename,nrows = maxAge+1)
  tmp$fullRVoptions<-(tmp$Client1==1& tmp$Client2==0) | 
    (tmp$Client1==0 & tmp$Client2==1)
  tmptimeInter<-
    tmp[fullRVoptions==TRUE,.(Prob.RV.V=mean(Choice)),
      by=.(Interv=floor(Age/interV),Training,Alpha,
           Gamma,Neta,Outbr,AlphaTh)]
  if(length(extPar)>0){
    tmptimeInter[,eval(extPar):=parVal]
  }
  return(tmptimeInter)
}

# Computes RV preference for a certain final proportion of the learning trial --

file2lastProp<-function(filename,prop,outPar=NULL,genfold=NULL,full.path=FALSE)
{
  if(length(outPar)>0){
    extPar<-grep(outPar,strsplit(filename,"_/")[[1]],
               value=TRUE)
    parVal<-as.numeric(gsub("[[:alpha:]]",extPar,replacement = ''))
    extPar<-gsub("[[:digit:]]",extPar,replacement = '')
  }
  if(is.null(genfold)) 
    if(full.path)   tmp<-fread(filename)
    else tmp<-fread(here("Simulations",filename))
    else tmp<-fread(here("Simulations",genfold,filename))
  # resPtmp<-as.numeric(gsub("[[:alpha:]]", "", grep('pR',strsplit(filename,'_')[[1]],value=TRUE)))*0.1
  # visPtmp<-as.numeric(gsub("[[:alpha:]]", "", grep('pV',strsplit(filename,'_')[[1]],value=TRUE)))*0.1
  # tmp$resProb<-rep(resPtmp,dim(tmp)[1])
  # tmp$visProb<-rep(visPtmp,dim(tmp)[1])
  tmp$fullRVoptions<-(tmp$Client1==1& tmp$Client2==0) | (tmp$Client1==0 & tmp$Client2==1)
  tmp<-tmp[Age>max(Age)*prop]
  lastchunk<-tmp[fullRVoptions==TRUE,
                    list(Prob.RV.V=mean(Choice)),
                    by=.(Training,Gamma,Neta,pR,pV,Outbr)]
  if(length(outPar)>0){
    lastchunk[,eval(outPar):=parVal]
  }
  return(lastchunk)
}

# Logistic function to compute probability from preference ---------------------

logist<-function(theta1,theta2){
  return(1/(1+exp(-(theta1-theta2))))
}

# Visualize difference in parameters between 2 JSON files

diffJsons<-function(json1,json2){
  print("JSON.1")
  print(unlist(json1)[unlist(json1)!=unlist(json2)])
  print("JSON.2")
  print(unlist(json2)[unlist(json1)!=unlist(json2)])
}

# Compute the number of rounds it takes to develop a certain preference fo V ---

loadDataFirstReach<-function(filename,bound){
  tmp<-fread(filename)
  tmp$fullRVoptions<-(tmp$Client1==1& tmp$Client2==0) | (tmp$Client1==0 & tmp$Client2==1)
  tmp[,ReachedCut:=logist(theta1 = ThetaV,theta2 = ThetaR)>bound]
  tmpFirstR<-tmp[ReachedCut==TRUE,.(firstReach=min(Age),
                                       Prob.RV.V=logist(
                                         theta1 = ThetaV[Age==min(Age)],
                                         theta2 = ThetaR[Age==min(Age)])),
                    by=.(Training,Gamma,Neta,pR,pV,Outbr)]
  return(tmpFirstR)
}


# function that creates new folders for simulations results --------------------
check_create.dir<-function(folder,param,values){
  setwd(folder)
  listfolders<-paste(param,values,"_",sep = "")  
  currFolders<-lapply(listfolders,dir.exists)
  if(sum(currFolders>0)){
    warning("At least one of the folders already exists \n Please check",
            immediate. = TRUE)
    print(cbind(listfolders,currFolders))
    # ans<-readline("Want to continue?")
    # if(substr(ans, 1, 1) == "y"){
    #   lapply(listfolders,dir.create)
    #   return(listfolders)
    # }
    # else{
    #   return(listfolders)
    # }
  }else{
    lapply(listfolders,dir.create)
    return(listfolders)
  }
}

loadMCMCrun<-function(scen,burn.in=1000,thinning=100){
  listFiles<-list.files(here("Simulations",scen))
  runs<-grep("MCMCchain",listFiles,value = TRUE)
  raw<-do.call(rbind,lapply(runs, function(file){
    rundata<-fread(here("Simulations",scen,file))
    seedNum<-strsplit(file,"_")[[1]][3]
    seedNum<-as.numeric(gsub("[[:alpha:]]",seedNum,replacement = ''))
    rundata[,seed:=rep(seedNum,dim(rundata)[1])]
  }))
  filtered<-raw[iteration>burn.in & iteration %% thinning==0,]
  return(filtered)
}
  
  
  
  