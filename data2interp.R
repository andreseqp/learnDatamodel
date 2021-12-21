######################### Abundance functions ##################################

# Function to create the interpolation -----------------------------------------

AbundData2interp<-function(dataAbun,
                           npoints=100, # size of the interpolations
                           Var2int      # variable to interpolate
                           ){

  # interpolate real data
  interpData<-with(dataAbun,
                      {interp(x=pR,y=pV,z=get(Var2int),duplicate = "mean",
                              nx=npoints,ny=npoints)})
  # Create structure to transpose interpolated data
  interpDataTrans<-data.table(matrix(0,nrow = npoints*npoints,ncol = 4))
  names(interpDataTrans)<-c("resProb","visProb",Var2int,"notProb")
  # Transpose interpolated data
  for (i in 1:npoints) {
    for (j in 1:npoints) {
      interpDataTrans[(i-1)*npoints+j,resProb:=interpData$x[i]]
      interpDataTrans[(i-1)*npoints+j,visProb:=interpData$y[j]]
      interpDataTrans[(i-1)*npoints+j,eval(Var2int):=interpData$z[i,j]]
    }
  }
  # Complementary probability - "absence"
  interpDataTrans[,4]<-1-interpDataTrans[,1]-interpDataTrans[,2]
  # Select relavant data
  interpDataTrans<-interpDataTrans[resProb+visProb<1]
  return(interpDataTrans)
}

# Function to interpolate abundance vs leaving probability --------------------- 

AbundLeavData2interp<-function(dataAbun,
                           npoints=100, # size of the interpolations
                           Var2int ){     # variable to interpolate

  
  # interpolate real data
  interpData<-with(dataAbun,
                   {interp(x=pA,y=Vlp,z=get(Var2int),duplicate = "mean",
                           nx=npoints,ny=npoints)})
  # Create structure to transpose interpolated data
  interpDataTrans<-data.table(matrix(0,nrow = npoints*npoints,ncol = 3))
  names(interpDataTrans)<-c("pAbs","VLeavProb",Var2int)
  # Transpose interpolated data
  for (i in 1:npoints) {
    for (j in 1:npoints) {
      interpDataTrans[(i-1)*npoints+j,pAbs:=interpData$x[i]]
      interpDataTrans[(i-1)*npoints+j,VLeavProb:=interpData$y[j]]
      interpDataTrans[(i-1)*npoints+j,eval(Var2int):=interpData$z[i,j]]
    }
  }
  return(interpDataTrans)
}


