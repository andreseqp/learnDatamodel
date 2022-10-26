library(Rcpp)
library(here)
library("BayesianTools")
library(data.table)
library("jsonlite")

sourceCpp("test.cpp")

timesTwo(3)

?VSEM

focal<-new(model_param)
focal$set_gamma(0.5,0.5)
focal$logist()
focal$methodJson()

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

sourceCpp("ActCrit_R.cpp")


foc.param<-list(alphaC=0.01,alphaA=0.01,scaleConst=150,
     gamma=c(0.1,0.2),negReward=c(0.5,0.8),
     probFAA=c(1,1),interpReg=1,slopRegRelAC=2,
     slopRegPVL=2)

fieldData<-fread(here("Data","data_cleaner_abs_threa1.5.txt"))
names(fieldData)[4:8]<-c("abund_clean","abund_visitors","abund_resid",
                         "prob_Vis_Leav","group")

scenario<-"testBayesianTools"

typeof(outParam)

param_mcmc$agentScen <-0
outParam<-toJSON(param_mcmc,auto_unbox = TRUE,pretty = TRUE)


paramAddr<-here("Simulations",paste0(scenario,"_"),
                "parametersMCMCcluster_2.json")

predicton<-do_simulation(emp_data = fieldData,
              focal_param = foc.param,
              fileStr= outParam)
# ,
#               fileStr =   paramAddr)
