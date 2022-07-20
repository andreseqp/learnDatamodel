##### Analysis of the Aproximate Bayesian Computation to estimating learning
##### parameters in the AC models using data from experiments and the field

## Libraries
library(data.table)
library(coda)
library(ggplot2)
library(cowplot)
library(ggdist)
library(ggpubr)
source(here("../R_files/posPlots.R"))
source(here("loadData.R"))
library("jsonlite")

# Name of the folders containing the MCMC runs
# scen1<-"MCMCclean_gam_nRew_sca_"
# scen2<-"MCMCclean_gam_sca_"
# scen3<-"MCMCclean_Nrew_sca_"

scen1<-"MCMCclean_gam_Nrew_sca_"
scen2<-"MCMCclean_gam_sca_"
scen3<-"MCMCclean_Nrew_sca_"


# NLables names for the three scenarios
labels.Scen<-c("both","gam","Nrew")
labelsPlot.Scen<-c("Full Model","Chaining","Penalty")

## Load files --------------------------------------------------------------

# Defaults for MCMC analysis
burn.in<-1000
thinning<-100


MCMCdata1<-loadMCMCrun(scen1,thinning = thinning,burn.in = burn.in)
MCMCdata2<-loadMCMCrun(scen2,thinning = thinning,burn.in = burn.in)
MCMCdata3<-loadMCMCrun(scen3,thinning = thinning,burn.in = burn.in)


# Create list with the data, convenient to use CODA for diagnostics
mcmcList.1<-mcmc.list(do.call(list,lapply(unique(MCMCdata1$seed), function(repli){
  rundata<-MCMCdata1[seed==repli]
  mcmcRun<-mcmc(rundata[,.(gamma,negReward,scaleConst)],thin = thinning)#
  return(mcmcRun)
})))

mcmcList.2<-mcmc.list(do.call(list,lapply(unique(MCMCdata2$seed), function(repli){
  rundata<-MCMCdata2[seed==repli]
  mcmcRun<-mcmc(rundata[,.(gamma,scaleConst)],thin = thinning)
  return(mcmcRun)
})))


mcmcList.3<-mcmc.list(do.call(list,lapply(unique(MCMCdata3$seed), function(repli){
  rundata<-MCMCdata3[seed==repli]
  mcmcRun<-mcmc(rundata[,.(negReward,scaleConst)],thin = thinning)#
  return(mcmcRun)
})))

## Density plots with gg -------------------------------------------------------

## Define the boundaries for credible intervals using High Density Interval

intervals<-c(0.66,0.95)

cuts.1<-lapply(MCMCdata1[,.(gamma,negReward,scaleConst)],
               function(x){
                 c(mode_hdci(x,.width = intervals[c(2,1)])$ymin,
                   mode_hdci(x,.width = intervals)$ymax)
               })


cuts.2<-lapply(MCMCdata2[,.(gamma,scaleConst)],
               function(x){
                 c(mode_hdci(x,.width = intervals[c(2,1)])$ymin,
                   mode_hdci(x,.width = intervals)$ymax)
               })

cuts.3<-lapply(MCMCdata3[,.(negReward,scaleConst)],
               function(x){
                 c(mode_hdci(x,.width = intervals[c(2,1)])$ymin,
                   mode_hdci(x,.width = intervals)$ymax)
               })


# get all the likelihoods from the posterior for all models
loglikehoods.all<-data.table(lglikelihood=c(MCMCdata1$fit,MCMCdata2$fit,
                                            MCMCdata3$fit),
                             model=c(rep(labels.Scen[1],length(MCMCdata1$fit)),
                                     rep(labels.Scen[2],length(MCMCdata2$fit)),
                                     rep(labels.Scen[3],length(MCMCdata3$fit))
                             ))

# Read field data to calculate null log-likelihood
fieldData<-fread(here("Data","data_cleaner_abs_threa1.5.txt"))

nullLikehood<-sum(dbinom(x=fieldData[,score_visitor],
                         size = 20,prob = 0.5,log = TRUE))

# calculate pseudo R^2
loglikehoods.all[,Rsqrd:=1-lglikelihood/nullLikehood]


panel.rsqr.all<-
  ggplot(data=loglikehoods.all[is.finite(lglikelihood)],
         aes(y=model,x=Rsqrd,fill=model))+
  stat_halfeye(alpha=0.5,point_size=3)+
  theme_classic()+xlim(-0.5,0.25)+
  scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3"),
                    name = "", labels = labelsPlot.Scen)+
  ylab("")+
  scale_y_discrete(labels=labelsPlot.Scen)+
  geom_vline(xintercept = 0,color="black",size=0.5)+
  theme(legend.position = 'none')
# c(.32, .65)
# ,legend.key.size = unit(0.22,'cm'),
# axis.text.y = element_blank(),
# legend.text = element_text(size=7))

gam.both.post<-ggplot(data=MCMCdata1,aes(x=gamma)) +
  stat_halfeye(aes(fill=stat(cut(x,breaks = cuts.1$gamma))),
               point_interval = mode_hdi, .width = c(.66, .95),
               point_size=3,show.legend=FALSE) +
  labs(title=labelsPlot.Scen[1],
       subtitle = expression(gamma),x="",y="")+
  theme_classic()+
  scale_fill_manual(values=c("gray85","skyblue","gray85"))

nrew.both.post<-ggplot(data=MCMCdata1,aes(x=negReward)) +
  stat_halfeye(aes(fill=stat(cut(x,breaks = cuts.1$negReward))),
               point_interval = mode_hdi,
               .width = c(.66, .95),point_size=3,
               show.legend=FALSE)+theme_classic()+
  scale_fill_manual(values=c("gray85","skyblue","gray85"))+
  labs(title="",subtitle = expression(eta),x="",y="")+
  xlim(-2,2)
gam.gam.post<-ggplot(data=MCMCdata2,aes(x=gamma)) +
  stat_halfeye(aes(fill=stat(cut(x,breaks = cuts.2$gamma))),
               point_interval = mode_hdi, .width = c(.66, .95),
               point_size=3,
               show.legend=FALSE) +
  scale_fill_manual(values=c("gray85","skyblue","gray85"))+
  labs(title=labelsPlot.Scen[2],
       subtitle = expression(gamma),x="",y="")+
  theme_classic()
nrew.nrew.post<-ggplot(data=MCMCdata3,aes(x=negReward)) +
  stat_halfeye(aes(fill=stat(cut(x,breaks = cuts.3$negReward))),
               point_interval = mode_hdi, .width = c(.66, .95),
               point_size=3,
               show.legend=FALSE) +
  scale_fill_manual(values=c("gray85","skyblue","gray85"))+
  labs(title=labelsPlot.Scen[3],
       subtitle = expression(eta),x="",y="")+
  theme_classic()+xlim(-2,2)
