##### Analysis of the Aproximate Bayesian Computation to estimating learning
##### parameters in the AC models using data from experiments and the field

## Libraries
library(data.table)
library(coda)
library(ggplot2)
library(cowplot)
library(ggdist)
library(ggpubr)
library(here)
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
scen4<-"MCMCclean_PAA2_gam_Nrew_sca_"
scen5<-"MCMCclean_PAA2_gam_sca_"
scen6<-"MCMCclean_PAA2_Nrew_sca_"


# NLables names for the three scenarios
labels.Scen<-c("both","gam","Nrew")
labelsPlot.Scen<-c("Full Model","Chaining","Penalty")
labels.agente<-c("FAA","PAA")
## Load files --------------------------------------------------------------

# Defaults for MCMC analysis
burn.in<-1000
thinning<-100


MCMCdata1<-loadMCMCrun(scen1,thinning = thinning,burn.in = burn.in)
MCMCdata2<-loadMCMCrun(scen2,thinning = thinning,burn.in = burn.in)
MCMCdata3<-loadMCMCrun(scen3,thinning = thinning,burn.in = burn.in)
MCMCdata4<-loadMCMCrun(scen4,thinning = thinning,burn.in = burn.in)
MCMCdata5<-loadMCMCrun(scen5,thinning = thinning,burn.in = burn.in)
MCMCdata6<-loadMCMCrun(scen6,thinning = thinning,burn.in = burn.in)


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

MCMCdata4<-MCMCdata4[iteration<66800]

mcmcList.4<-mcmc.list(do.call(list,lapply(unique(MCMCdata4$seed), function(repli){
  rundata<-MCMCdata4[seed==repli]
  mcmcRun<-mcmc(rundata[,.(gamma,negReward,scaleConst)],thin = thinning)#
  return(mcmcRun)
})))

MCMCdata5<-MCMCdata5[iteration<MCMCdata5[,max(iteration),by=seed][,min(V1)]]

mcmcList.5<-mcmc.list(do.call(list,lapply(unique(MCMCdata5$seed), function(repli){
  rundata<-MCMCdata5[seed==repli]
  mcmcRun<-mcmc(rundata[,.(gamma,scaleConst)],thin = thinning)
  return(mcmcRun)
})))

MCMCdata6<-MCMCdata6[iteration<MCMCdata6[,max(iteration),by=seed][,min(V1)]]
mcmcList.6<-mcmc.list(do.call(list,lapply(unique(MCMCdata6$seed), function(repli){
  rundata<-MCMCdata6[seed==repli]
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

cuts.4<-lapply(MCMCdata4[,.(gamma,negReward,scaleConst)],
               function(x){
                 c(mode_hdci(x,.width = intervals[c(2,1)])$ymin,
                   mode_hdci(x,.width = intervals)$ymax)
               })


cuts.5<-lapply(MCMCdata5[,.(gamma,scaleConst)],
               function(x){
                 c(mode_hdci(x,.width = intervals[c(2,1)])$ymin,
                   mode_hdci(x,.width = intervals)$ymax)
               })

cuts.6<-lapply(MCMCdata6[,.(negReward,scaleConst)],
               function(x){
                 c(mode_hdci(x,.width = intervals[c(2,1)])$ymin,
                   mode_hdci(x,.width = intervals)$ymax)
               })


# get all the likelihoods from the posterior for all models
loglikehoods.all<-data.table(lglikelihood=c(MCMCdata1$fit,MCMCdata2$fit,
                                            MCMCdata3$fit
                                            ,MCMCdata4$fit,
                                            MCMCdata5$fit,MCMCdata6$fit
                                            ),
                             model=c(rep(labels.Scen[1],length(MCMCdata1$fit)),
                                     rep(labels.Scen[2],length(MCMCdata2$fit)),
                                     rep(labels.Scen[3],length(MCMCdata3$fit))
                                     ,rep(labels.Scen[1],length(MCMCdata4$fit)),
                                     rep(labels.Scen[2],length(MCMCdata5$fit)),
                                     rep(labels.Scen[3],length(MCMCdata6$fit))
                                     )
                             ,agent=rep(labels.agente,each=length(MCMCdata3$fit)*3)
                             )
loglikehoods.all[,model_agent:=
                   factor(paste(model,agent,sep = "_"),
                          levels=c("both_FAA","gam_FAA","Nrew_FAA",
                                   "both_PAA","gam_PAA","Nrew_PAA"),
                          labels = c("FAA full model", "FAA chaining",
                                     "FAA penalty","PAA full model",
                                     "PAA chaining","PAA penalty"))]

loglikehoods.all[,levels(model_agent)]

# Read field data to calculate null log-likelihood
fieldData<-fread(here("Data","data_cleaner_abs_threa1.5.txt"))

nullLikehood<-sum(dbinom(x=fieldData[,score_visitor],
                         size = 20,prob = 0.5,log = TRUE))

# Calculate pseudo R^2
loglikehoods.all[,Rsqrd:=1-lglikelihood/nullLikehood]


panel.rsqr.all<-
  ggplot(data=loglikehoods.all[is.finite(lglikelihood)],
         aes(y=model_agent,x=Rsqrd,fill=model))+
  stat_halfeye(alpha=0.5,point_size=3)+
  # facet_grid(~agent)+
  theme_classic()+xlim(-0.5,0.3)+
  scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3"),
                    name = "")+
  ylab("")+
  geom_vline(xintercept = 0,color="black",size=0.5)+
  theme(legend.position = 'none')
# c(.32, .65)
# ,legend.key.size = unit(0.22,'cm'),
# axis.text.y = element_blank(),
# legend.text = element_text(size=7))



# Posterior of gamma for all models

allGammapost<-data.table(
  gamma = c(MCMCdata1$gamma,MCMCdata2$gamma,MCMCdata4$gamma,MCMCdata5$gamma),
  model = factor(c(rep("FAA.both",length(MCMCdata1$gamma)),
                   rep("FAA.gamma",length(MCMCdata2$gamma)),
                   rep("PAA.both",length(MCMCdata4$gamma)),
                   rep("PAA.gamma",length(MCMCdata5$gamma))),
             labels = c("FAA full model", "FAA chaining",
                        "PAA full model","PAA chaining")),
  agent = factor(c(rep("FAA",length(MCMCdata1$gamma)+length(MCMCdata2$gamma)),
                   rep("PAA",length(MCMCdata4$gamma)+length(MCMCdata5$gamma))))
)

allGammapost.plot<-ggplot(allGammapost,aes(y=model,x=gamma))+
  stat_halfeye(aes(fill=agent,fill_ramp=stat(level)),.width = c(.66, .95),
               point_interval = mode_hdi, 
               point_size=1.5,stroke=3,position = "dodge") +
  labs(title=expression(gamma),fill_ramp="Interval",x="",y="")+
  theme_classic()+
  scale_fill_manual(values=c("dodgerblue2", "#E31A1C"))+
  scale_fill_ramp_discrete(na.translate = FALSE) +
  geom_vline(xintercept = 0)

alletapost<-data.table(
  negReward = c(MCMCdata1$negReward,MCMCdata3$negReward,
                MCMCdata4$negReward,MCMCdata6$negReward),
  model = factor(c(rep("FAA.both",length(MCMCdata1$negReward)),
                   rep("FAA.negReward",length(MCMCdata3$negReward)),
                   rep("PAA.both",length(MCMCdata4$negReward)),
                   rep("PAA.negReward",length(MCMCdata6$negReward))),
                 labels = c("FAA full model", "FAA penalty",
                            "PAA full model","PAA penalty")),
  agent = factor(c(rep("FAA",length(MCMCdata1$gamma)+length(MCMCdata2$gamma)),
                   rep("PAA",length(MCMCdata4$gamma)+length(MCMCdata5$gamma))))
)

alletapost.plot<-ggplot(alletapost,aes(y=model,x=negReward))+
  stat_halfeye(aes(fill=agent,fill_ramp=stat(level)),.width = c(.66, .95),
               point_interval = mode_hdi, 
               point_size=1.5,stroke=3,position = "dodgejust") +
  labs(title=expression(eta),fill_ramp="Interval",x="",y="")+
  theme_classic()+
  scale_fill_manual(values=c("dodgerblue2", "#E31A1C"))+
  scale_fill_ramp_discrete(na.translate = FALSE) +
  geom_vline(xintercept = 0)+
  xlim(-1,13)
  



gam.both.post.FAA<-ggplot(data=MCMCdata1,aes(x=gamma)) +
  stat_halfeye(aes(fill=stat(level)),#cut(x,breaks = cuts.1$gamma))),
               point_interval = mode_hdi, .width = c(.66, .95),
               point_size=3,show.legend=FALSE) +
  labs(title=labelsPlot.Scen[1],
       subtitle = expression(gamma),x="",y="")+
  theme_classic()+
  scale_fill_manual(values=c("gray85","skyblue","gray85"))

gam.both.post.PAA<-ggplot(data=MCMCdata4,aes(x=gamma)) +
  stat_halfeye(aes(fill=stat(cut(x,breaks = cuts.4$gamma))),
               point_interval = mode_hdi, .width = c(.66, .95),
               point_size=3,show.legend=FALSE) +
  labs(title=labelsPlot.Scen[1],
       subtitle = expression(gamma),x="",y="")+
  theme_classic()+
  scale_fill_manual(values=c("gray85","skyblue","gray85"))


nrew.both.post.FAA<-ggplot(data=MCMCdata1,aes(x=negReward)) +
  stat_halfeye(aes(fill=stat(cut(x,breaks = cuts.1$negReward))),
               point_interval = mode_hdi,
               .width = c(.66, .95),point_size=3,
               show.legend=FALSE)+theme_classic()+
  scale_fill_manual(values=c("gray85","skyblue","gray85"))+
  labs(title="",subtitle = expression(eta),x="",y="")+
  xlim(-2,2)

nrew.both.post.PAA<-ggplot(data=MCMCdata4,aes(x=negReward)) +
  stat_halfeye(aes(fill=stat(cut(x,breaks = cuts.5$negReward))),
               point_interval = mode_hdi,
               .width = c(.66, .95),point_size=3,
               show.legend=FALSE)+theme_classic()+
  scale_fill_manual(values=c("gray85","skyblue","gray85"))+
  labs(title="",subtitle = expression(eta),x="",y="")+
  xlim(-2,2)


gam.gam.post.FAA<-ggplot(data=MCMCdata2,aes(x=gamma)) +
  stat_halfeye(aes(fill=stat(cut(x,breaks = cuts.2$gamma))),
               point_interval = mode_hdi, .width = c(.66, .95),
               point_size=3,
               show.legend=FALSE) +
  scale_fill_manual(values=c("gray85","skyblue","gray85"))+
  labs(title=labelsPlot.Scen[2],
       subtitle = expression(gamma),x="",y="")+
  theme_classic()

gam.gam.post.PAA<-ggplot(data=MCMCdata5,aes(x=gamma)) +
  stat_halfeye(aes(fill=stat(cut(x,breaks = cuts.5$gamma))),
               point_interval = mode_hdi, .width = c(.66, .95),
               point_size=3,
               show.legend=FALSE) +
  scale_fill_manual(values=c("gray85","skyblue","gray85"))+
  labs(title=labelsPlot.Scen[2],
       subtitle = expression(gamma),x="",y="")+
  theme_classic()


nrew.nrew.post.FAA<-ggplot(data=MCMCdata3,aes(x=negReward)) +
  stat_halfeye(aes(fill=stat(cut(x,breaks = cuts.3$negReward))),
               point_interval = mode_hdi, .width = c(.66, .95),
               point_size=3,
               show.legend=FALSE) +
  scale_fill_manual(values=c("gray85","skyblue","gray85"))+
  labs(title=labelsPlot.Scen[3],
       subtitle = expression(eta),x="",y="")+
  theme_classic()+xlim(-2,2)

nrew.nrew.post.PAA<-ggplot(data=MCMCdata6,aes(x=negReward)) +
  stat_halfeye(aes(fill=stat(cut(x,breaks = cuts.6$negReward))),
               point_interval = mode_hdi, .width = c(.66, .95),
               point_size=3,
               show.legend=FALSE) +
  scale_fill_manual(values=c("gray85","skyblue","gray85"))+
  labs(title=labelsPlot.Scen[3],
       subtitle = expression(eta),x="",y="")+
  theme_classic()+xlim(-2,2)
