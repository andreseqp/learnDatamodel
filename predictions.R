#### Figures showing the fit between model and data ############################
#----------- Figure 2 
##
library(data.table)
library(ggplot2)
library("RColorBrewer")
library("patchwork") 
library(ggpubr)
library(ggdist)
library(ggforce)
library(concaveman)
library(ggalt)
library(here)
library(cowplot)
source(here("loadData.R"))
source(here("data2interp.R"))
require('akima')
source(here("aesth_par.R"))
source(here("../R_files/posPlots.R"))


# Folder where the three models are
# scen1<-"MCMCclean_gam_Nrew_sca_"
# scen2<-"MCMCclean_gam_sca_"
# scen3<-"MCMCclean_Nrew_sca_"

# Get file with the predictions using the mode of the posterior
predfileMode.both<-grep(".txt",grep("round",list.files(here("Simulations",
                    scen1)),value = T),value = T)
predfileMode.gam<-grep(".txt",grep("round",list.files(here("Simulations",
                    scen2)),value = T),value = T)
predfileMode.Nrew<-grep(".txt",grep("round",list.files(here("Simulations",
                    scen3)),value = T),value = T)

# Get predictions from sampling the posterior 
predfileSamples.both<-grep(".txt",grep("round",list.files(here("Simulations",
                                  scen1,"samplesPost_"),
                                  full.names = TRUE),value = T), value = T)
predfileSamples.gam<-grep(".txt",grep("round",list.files(here("Simulations",
                                scen2,"samplesPost_"),
                                full.names = TRUE),value = T),value = T)
predfileSamples.Nrew<-grep(".txt",grep("round",list.files(here("Simulations",
                                scen3,"samplesPost_"),
                               full.names = TRUE),value = T),value = T)

# Load predictions with the mode of the posterior
predictDataMode.both<-fread(here("Simulations",scen1,
                      predfileMode.both))
predictDataMode.gam<-fread(here("Simulations",scen2,
                             predfileMode.gam))
predictDataMode.Nrew<-fread(here("Simulations",scen3,
                                predfileMode.Nrew))

# Load predictions with the samples from the posterior
predictDataSamps.both<-do.call(rbind,lapply(predfileSamples.both,fread))
predictDataSamps.both[,id_samp:=rep(1:length(predfileSamples.both),each=120)]

predictDataSamps.gam<-do.call(rbind,lapply(predfileSamples.gam,fread))
predictDataSamps.gam[,id_samp:=rep(1:length(predfileSamples.gam),each=120)]

predictDataSamps.Nrew<-do.call(rbind,lapply(predfileSamples.Nrew,fread))
predictDataSamps.Nrew[,id_samp:=rep(1:length(predfileSamples.Nrew),each=120)]



# SUmmarize data sets by location
fieldatabyLoc.both<-predictDataMode.both[,.(probVisi.data=mean(visitorChoices)/20,
                              probvisitor.pred=max(visitorChoices_pred),
                              re.abund.clean=max(rel.abund.cleaners),
                              prob.Vis.leave=max(prob.Vis.Leav),
                              model="both"),by=site_year]

fieldatabyLoc.gam<-predictDataMode.gam[,.(probVisi.data=mean(visitorChoices)/20,
                                   probvisitor.pred=max(visitorChoices_pred),
                                   re.abund.clean=max(rel.abund.cleaners),
                                   prob.Vis.leave=max(prob.Vis.Leav),
                                   model="gamma"
                                   ),by=site_year]

fieldatabyLoc.Nrew<-predictDataMode.Nrew[,.(probVisi.data=mean(visitorChoices)/20,
                                          probvisitor.pred=max(visitorChoices_pred),
                                          re.abund.clean=max(rel.abund.cleaners),
                                          prob.Vis.leave=max(prob.Vis.Leav),
                                          model="eta"
                                          ),by=site_year]

fieldatabyLocSamps.both<-predictDataSamps.both[,.(probVisi.data=mean(visitorChoices)/20,
                                        probvisitor.pred=max(visitorChoices_pred),
                                        re.abund.clean=max(rel.abund.cleaners),
                                        prob.Vis.leave=max(prob.Vis.Leav)),
                                        by=.(site_year,id_samp)]

fieldatabyLocSamps.gam<-predictDataSamps.gam[,.(probVisi.data=mean(visitorChoices)/20,
                                      probvisitor.pred=max(visitorChoices_pred),
                                      re.abund.clean=max(rel.abund.cleaners),
                                      prob.Vis.leave=max(prob.Vis.Leav)),
                                      by=.(site_year,id_samp)]

fieldatabyLocSamps.Nrew<-predictDataSamps.Nrew[,.(probVisi.data=mean(visitorChoices)/20,
                                                probvisitor.pred=max(visitorChoices_pred),
                                                re.abund.clean=max(rel.abund.cleaners),
                                                prob.Vis.leave=max(prob.Vis.Leav)),
                                             by=.(site_year,id_samp)]

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0,1))


## Load data for the contour plot using the mode of the posterior --------------

pointInt<-200

# Full model
simsDir.both<-here("Simulations",scen1)

FIAlastQuarData<-do.call(rbind,lapply(
  getFilelist(simsDir.both,fullNam = TRUE)$p1,file2lastProp,0.70,outPar="Vlp",
  full.path=TRUE))


FIAlastQuarData[,pA:=1-pR-pV]

FIA.stats<-FIAlastQuarData[,.(meanProb=mean(Prob.RV.V),
                              upIQR=fivenum(Prob.RV.V)[4],
                              lowIQR=fivenum(Prob.RV.V)[2])
                           ,by=.(Neta,Gamma,pR,pV,Vlp)]

FIA.stats$pA<-round(1-FIA.stats$pR-FIA.stats$pV,1)

FIAinterpData.both<-AbundLeavData2interp(FIAlastQuarData,
                                    Var2int = "Prob.RV.V",npoints = pointInt)

# Only Gamma

simsDir.gam<-here("Simulations",scen2)

FIAlastQuarData<-do.call(rbind,lapply(
  getFilelist(simsDir.gam,fullNam = TRUE)$p1,file2lastProp,0.70,outPar="Vlp",
  full.path=TRUE))

FIAlastQuarData[,pA:=1-pR-pV]

FIA.stats<-FIAlastQuarData[,.(meanProb=mean(Prob.RV.V),
                              upIQR=fivenum(Prob.RV.V)[4],
                              lowIQR=fivenum(Prob.RV.V)[2])
                           ,by=.(Neta,Gamma,pR,pV,Vlp)]

FIA.stats$pA<-round(1-FIA.stats$pR-FIA.stats$pV,1)

FIAinterpData.gam<-AbundLeavData2interp(FIAlastQuarData,
                                    Var2int = "Prob.RV.V",npoints = pointInt)

# Only eta

simsDir.Nrew<-here("Simulations",scen3)

FIAlastQuarData<-do.call(rbind,lapply(
  getFilelist(simsDir.Nrew,fullNam = TRUE)$p1,file2lastProp,0.70,outPar="Vlp",
  full.path=TRUE))

FIAlastQuarData[,pA:=1-pR-pV]

FIA.stats<-FIAlastQuarData[,.(meanProb=mean(Prob.RV.V),
                              upIQR=fivenum(Prob.RV.V)[4],
                              lowIQR=fivenum(Prob.RV.V)[2])
                           ,by=.(Neta,Gamma,pR,pV,Vlp)]

FIA.stats$pA<-round(1-FIA.stats$pR-FIA.stats$pV,1)

FIAinterpData.Nrew<-AbundLeavData2interp(FIAlastQuarData,
                                        Var2int = "Prob.RV.V",npoints = pointInt)

# rm(list("FIAlastQuarData","FIA.stats"))


## Finally!! plot the predictions ----------------------------------------------

# For convenience change names 
names(FIAinterpData.both)<-c("rel.abund.cleaners","prob.Vis.Leav","market_binomial_data")
names(fieldatabyLoc.both)[c(4,5,2)]<-c("rel.abund.cleaners","prob.Vis.Leav","market_binomial_data")

names(FIAinterpData.gam)<-c("rel.abund.cleaners","prob.Vis.Leav","market_binomial_data")
names(fieldatabyLoc.gam)[c(4,5,2)]<-c("rel.abund.cleaners","prob.Vis.Leav","market_binomial_data")

names(FIAinterpData.Nrew)<-c("rel.abund.cleaners","prob.Vis.Leav","market_binomial_data")
names(fieldatabyLoc.Nrew)[c(4,5,2)]<-c("rel.abund.cleaners","prob.Vis.Leav","market_binomial_data")


# Calculate Log-likelihoods
predictDataMode.both[,log.like:=dbinom(visitorChoices,
                      size = 20,prob = visitorChoices_pred,log = TRUE)]

predictDataMode.gam[,log.like:=dbinom(visitorChoices,
                     size = 20,prob = visitorChoices_pred,log = TRUE)]

predictDataMode.Nrew[,log.like:=dbinom(visitorChoices,
                    size = 20,prob = visitorChoices_pred,log = TRUE)]

# Calculate pseudo R^2
rsqr.both.McFadden<-1-predictDataMode.both[,sum(log.like)]/
  sum(dbinom(x=predictDataMode.both[,visitorChoices],size = 20,prob = 0.5,log = TRUE))

rsqr.gam.McFadden<-1-predictDataMode.gam[,sum(log.like)]/
  sum(dbinom(x=predictDataMode.gam[,visitorChoices],size = 20,prob = 0.5,log = TRUE))

rsqr.Nrew.McFadden<-1-predictDataMode.Nrew[,sum(log.like)]/
  sum(dbinom(x=predictDataMode.Nrew[,visitorChoices],size = 20,prob = 0.5,log = TRUE))

axislabSize<-8
axisSize<-12
margPlots<-unit(c(0,0,0,0),"cm")
# Panel full model contour
cont.obs.pred.both<- ggplot(data = FIAinterpData.both,aes(x=rel.abund.cleaners,y=prob.Vis.Leav,
                                                fill=market_binomial_data))+
 geom_raster(interpolate = TRUE) +  
  scale_fill_gradientn(limits=c(0.3,1),colours= myPalette(100))+theme_classic()+
  geom_point(data = fieldatabyLoc.both,aes(fill=market_binomial_data),size=5,
             shape=21,color="black")+sc+xlab("")+
  ylab("Prob. of visitor leaving")+
  labs(fill="Probability \n of choosing a visitor")+#
  theme(legend.position ="bottom",
        axis.text = element_text(size=axisSize),
        axis.title.x = element_text(size=axislabSize),
        axis.title.y = element_text(size=axislabSize),
        plot.margin = margPlots)

# Panel full model scatter
# 
scatter.obs.pred.both<-ggplot(data = fieldatabyLoc.both,
                              aes(y=market_binomial_data,x=probvisitor.pred,
                                  colour=site_year))+
  # stat_gradientinterval(aes(slab_alpha=F),point_interval = mode_hdi,point_size=2,
  #                       fill_type = "segments")+
  geom_abline(slope=1)+ylab("Observed")+xlab("")+
  geom_point(shape=15,size=2)+
  geom_mark_ellipse(aes(color=model))+
  ylim(0.4,0.9)+xlim(0.45,0.85)+
  guides(color=guide_legend(title=""))+
  scale_color_manual(values = c("black",multDiscrPallet[1:12]))+
  theme_classic()+theme(legend.position ="top",#c(.9, .4),
        axis.text = element_text(size=axisSize),
        axis.title.x = element_text(size=axislabSize),
        axis.title.y = element_text(size=axislabSize),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.1,'cm'),
        axis.text.x = element_blank(),
        plot.margin = margPlots)+
  guides(colour=guide_legend(ncol=3,title="",override.aes = list(size=1)))

  
  

# Panel future reward contour
cont.obs.pred.gam<- ggplot(data = FIAinterpData.gam,aes(x=rel.abund.cleaners,y=prob.Vis.Leav,
                                                          fill=market_binomial_data))+
  geom_raster(interpolate = TRUE) +  
  scale_fill_gradientn(limits=c(0.3,1),colours= myPalette(100))+theme_classic()+
  geom_point(data = fieldatabyLoc.gam,aes(fill=market_binomial_data),size=5,
             shape=21,color="black")+sc+xlab("Relative cleaner abundance")+
  ylab("Prob. of visitor leaving")+
  labs(fill="Probability \n of choosing \n a visitor")+
  theme(legend.position ="none",
    axis.text = element_text(size=axisSize),
    axis.title.x = element_text(size=axislabSize),
    axis.title.y = element_text(size=axislabSize),
    plot.margin =  margPlots)

# Panel future reward scatter
# ggplot(data=fieldatabyLocSamps.gam,
#                              aes(y=probVisi.data,colour=site_year,x=probvisitor.pred))+
  # stat_gradientinterval(aes(slab_alpha=F),point_interval = mode_hdi,point_size=2,
                        # fill_type = "segments")+
scatter.obs.pred.gam<-ggplot(data=fieldatabyLoc.gam,aes(y=market_binomial_data,colour=site_year,
                                  x=probvisitor.pred))+
  geom_abline(slope=1)+ylab("Observed")+xlab("")+
  geom_point(shape=15,size=2)+
  geom_mark_ellipse(aes(color=model))+
  ylim(0.4,0.9)+xlim(0.45,0.85)+
  guides(color=guide_legend(title="Location"))+
  scale_color_manual(values = c("black",multDiscrPallet[1:12]))+
  theme_classic()+
  theme(legend.position ="none",#c(.85, .4),#
        legend.key.size = unit(0.20,'cm'),
        axis.text = element_text(size=axisSize),
        axis.title.x = element_text(size=axislabSize),
        axis.title.y = element_text(size=axislabSize),
        axis.text.x = element_blank(),
        plot.margin = margPlots)+
  guides(colour=guide_legend(ncol=1))

# Panel negative reward contour
cont.obs.pred.Nrew<- ggplot(data = FIAinterpData.Nrew,aes(x=rel.abund.cleaners,y=prob.Vis.Leav,
                                                        fill=market_binomial_data))+
  geom_raster(interpolate = TRUE) +  
  scale_fill_gradientn(limits=c(0.3,1),colours= myPalette(100))+theme_classic()+
  geom_point(data = fieldatabyLoc.gam,aes(fill=market_binomial_data),size=5,
             shape=21,color="black")+sc+xlab("Relative cleaner abundance")+
  ylab("Prob. of visitor leaving")+
  labs(fill="Probability \n of choosing \n a visitor")+
  theme(legend.position = "none",
    axis.text = element_text(size=axisSize),
    axis.title.x = element_text(size=axislabSize),
    axis.title.y = element_text(size=axislabSize),
    plot.margin = margPlots)

# Panel negative reward scatter


# ggplot(data=fieldatabyLocSamps.Nrew,
#                               aes(y=probVisi.data,colour=site_year,x=probvisitor.pred))+
  # stat_gradientinterval(aes(slab_alpha=F),point_interval = mode_hdi,point_size=2,
  #                       fill_type = "segments")+
scatter.obs.pred.Nrew<-
  ggplot(data = fieldatabyLoc.Nrew,aes(y=market_binomial_data,x=probvisitor.pred,
                                     color=site_year))+
  geom_abline(slope=1)+ylab("Observed")+xlab("Predicted")+
  geom_point(shape=15,size=2)+
  geom_mark_ellipse(aes(color=model))+
  ylim(0.4,0.9)+xlim(0.45,0.85)+
  guides(color=guide_legend(title="Location"))+
  scale_color_manual(values = c("black",multDiscrPallet[1:12]),
                     name = "Location")+
  theme_classic()+theme(legend.position ="none",#c(.9, .4),
                      axis.text = element_text(size=axisSize),
                      axis.title.x = element_text(size=axislabSize),
                      axis.title.y = element_text(size=axislabSize),
                      legend.text = element_text(size=6),
                      legend.key.size = unit(0.1,'cm'),
                      plot.margin = margPlots)+
  guides(colour=guide_legend(ncol=3,title="",override.aes = list(size=1)))
  
ScatterData.All<-rbind(fieldatabyLoc.both,fieldatabyLoc.gam,fieldatabyLoc.Nrew)
ScatterData.All[,model:=factor(model,levels = c("both","gamma","eta"),
                               labels = c("Full model",
                                          "Chaining","Penalty"))]


scatter.all.plot<-
  ggplot(data = ScatterData.All,aes(y=market_binomial_data,
                                    x=probvisitor.pred,
                                       color=model))+
  geom_abline(slope=1)+ylab("Observed")+xlab("Predicted")+
  geom_point(shape=15,size=2)+
  geom_encircle()+
  ylim(0.4,0.9)+xlim(0.45,0.85)+
  scale_color_manual(values = c("black",multDiscrPallet[1:12]))+
  theme_classic()+theme(legend.position ="top",
                        axis.text = element_text(size=axisSize),
                        axis.title.x = element_text(size=axislabSize),
                        axis.title.y = element_text(size=axislabSize),
                        legend.text = element_text(size=10),
                        legend.key.size = unit(0.8,'cm'),
                        plot.margin = margPlots)+
  guides(colour=guide_legend(ncol=3,title="Model"),override.aes = list(size=2))


