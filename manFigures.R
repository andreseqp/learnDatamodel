#### Figures showing the fit between model and data ############################
#----------- Figure 2 
##
library(data.table)
library(ggplot2)
library("RColorBrewer")
library("patchwork") 
library(ggpubr)
library(ggdist)
library(here)
library(cowplot)
source(here("loadData.R"))
source(here("data2interp.R"))
require('akima')
source(here("aesth_par.R"))
source(here("../R_files/posPlots.R"))

# Folder where the three models are
scenario<-"MCMCclean_gam_Nrew_sca"
scenario2<-"MCMCclean_gam_sca"
scenario3<-"MCMCclean_Nrew_sca"

# Get file with the predictions using the mode of the posterior
predfileMode.both<-grep(".txt",grep("round",list.files(here("Simulations",paste0(scenario,"_"))),value = T),
                    value = T)
predfileMode.gam<-grep(".txt",grep("round",list.files(here("Simulations",paste0(scenario2,"_"))),value = T),
                   value = T)
predfileMode.Nrew<-grep(".txt",grep("round",list.files(here("Simulations",paste0(scenario3,"_"))),value = T),
                       value = T)

# Get predictions from sampling the posterior 
predfileSamples.both<-grep(".txt",grep("round",list.files(here("Simulations",
                                  paste0(scenario,"_"),"samplesPost_"),
                                  full.names = TRUE),value = T), value = T)
predfileSamples.gam<-grep(".txt",grep("round",list.files(here("Simulations",
                                paste0(scenario2,"_"),"samplesPost_"),
                                full.names = TRUE),value = T),value = T)
predfileSamples.Nrew<-grep(".txt",grep("round",list.files(here("Simulations",
                                paste0(scenario3,"_"),"samplesPost_"),
                               full.names = TRUE),value = T),value = T)

# Load predictions with the mode of the posterior
predictDataMode.both<-fread(here("Simulations",paste0(scenario,"_"),
                      predfileMode.both))
predictDataMode.gam<-fread(here("Simulations",paste0(scenario2,"_"),
                             predfileMode.gam))
predictDataMode.Nrew<-fread(here("Simulations",paste0(scenario3,"_"),
                                predfileMode.Nrew))

# Load predictions with the samples from the posterior
predictDataSamps.both<-do.call(rbind,lapply(predfileSamples.both,fread))
predictDataSamps.both[,id_samp:=rep(1:100,each=120)]

predictDataSamps.gam<-do.call(rbind,lapply(predfileSamples.gam,fread))
predictDataSamps.gam[,id_samp:=rep(1:100,each=120)]

predictDataSamps.Nrew<-do.call(rbind,lapply(predfileSamples.gam,fread))
predictDataSamps.Nrew[,id_samp:=rep(1:100,each=120)]



# SUmmarize data sets by location
fieldatabyLoc.both<-predictDataMode.both[,.(probVisi.data=mean(visitorChoices)/20,
                              probvisitor.pred=max(visitorChoices_pred),
                              re.abund.clean=max(rel.abund.cleaners),
                            prob.Vis.leave=max(prob.Vis.Leav)),by=site_year]

fieldatabyLoc.gam<-predictDataMode.gam[,.(probVisi.data=mean(visitorChoices)/20,
                                   probvisitor.pred=max(visitorChoices_pred),
                                   re.abund.clean=max(rel.abund.cleaners),
                                   prob.Vis.leave=max(prob.Vis.Leav)),by=site_year]

fieldatabyLoc.Nrew<-predictDataMode.Nrew[,.(probVisi.data=mean(visitorChoices)/20,
                                          probvisitor.pred=max(visitorChoices_pred),
                                          re.abund.clean=max(rel.abund.cleaners),
                                          prob.Vis.leave=max(prob.Vis.Leav)),by=site_year]

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
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0,0.8))


## Load data for the contour plot using the mode of the posterior --------------

# Full model
simsDir.both<-here("Simulations",paste0(scenario,"_"))

FIAlastQuarData<-do.call(rbind,lapply(
  getFilelist(simsDir.both,fullNam = TRUE)$p1,file2lastProp,0.70,outPar="Vlp",
  full.path=TRUE))


FIAlastQuarData[,pA:=1-pR-pV]

FIA.stats<-FIAlastQuarData[,.(meanProb=mean(Prob.RV.V),
                              upIQR=fivenum(Prob.RV.V)[4],
                              lowIQR=fivenum(Prob.RV.V)[2])
                           ,by=.(Neta,Gamma,pR,pV,Vlp)]

FIA.stats$pA<-round(1-FIA.stats$pR-FIA.stats$pV,1)

pointInt<-200

FIAinterpData.both<-AbundLeavData2interp(FIAlastQuarData,
                                    Var2int = "Prob.RV.V",npoints = pointInt)

# Only Gamma

simsDir.gam<-here("Simulations",paste0(scenario2,"_"))

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

simsDir.Nrew<-here("Simulations",paste0(scenario3,"_"))

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
# Panel full model contour
cont.obs.pred.both<- ggplot(data = FIAinterpData.both,aes(x=rel.abund.cleaners,y=prob.Vis.Leav,
                                                fill=market_binomial_data))+
 geom_raster(interpolate = TRUE) +  
  scale_fill_gradientn(limits=c(0.3,1),colours= myPalette(100))+theme_classic()+
  geom_point(data = fieldatabyLoc.both,aes(fill=market_binomial_data),size=5,
             shape=21,color="black")+sc+xlab("")+
  ylab("Probability of visitor leaving")+
  labs(fill="")+#Probability \n of choosing \n a visitor
  theme(legend.position ="top",
        axis.text = element_text(size=axisSize),
        axis.title.x = element_text(size=axislabSize),
        axis.title.y = element_text(size=axislabSize))

# Panel full model scatter
scatter.obs.pred.both<-ggplot(data=fieldatabyLocSamps.both,
                              aes(y=probVisi.data,colour=site_year,x=probvisitor.pred))+
  stat_gradientinterval(aes(slab_alpha=F),point_interval = mode_hdi,point_size=2)+
  geom_abline(slope=1)+ylab("Observed")+xlab("")+
  geom_point(data = fieldatabyLoc.both,aes(y=market_binomial_data,x=probvisitor.pred),
             shape=15,size=2)+
  ylim(0.45,0.85)+xlim(0.45,0.85)+
  guides(color=guide_legend(title=""))+
  scale_color_manual(values = multDiscrPallet[1:12])+
  theme_classic()+theme(legend.position ="top",#c(.9, .4),
        axis.text = element_text(size=axisSize),
        axis.title.x = element_text(size=axislabSize),
        axis.title.y = element_text(size=axislabSize),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.1,'cm'))+
  guides(colour=guide_legend(ncol=3,title="",override.aes = list(size=1)))

  
  # ggplot(data = fieldatabyLoc.both,aes(y=market_binomial_data,x=probvisitor.pred))+
  # geom_point(data = fieldatabyLocSamps.both,aes(x=probvisitor.pred,y=probVisi.data),
  #            color='grey')+
  # geom_point(size=4)+ylim(0.4,0.9)+xlim(0.4,0.9)+
  # geom_abline(slope=1)+ylab("Observed")+xlab("Predicted")+
  # ggtitle("")+
  # theme_classic()+
  # theme(plot.title = element_text(hjust = 0.5),
  #     axis.title.x = element_text(size=axisSize),axis.title.y = element_text(size=axislabSize),
  #     axis.text = element_text(size=14))+
  # geom_text(x = 0.8, y = 0.85, label = expression(y==x), parse = TRUE,size=3)+
  # geom_text(x = 0.8, y = 0.5, label = deparse(bquote(R^2==.(round(rsqr.both.McFadden,4)))), 
  #           parse = TRUE,size=3)

# Panel future reward contour
cont.obs.pred.gam<- ggplot(data = FIAinterpData.gam,aes(x=rel.abund.cleaners,y=prob.Vis.Leav,
                                                          fill=market_binomial_data))+
  geom_raster(interpolate = TRUE) +  
  scale_fill_gradientn(limits=c(0.3,1),colours= myPalette(100))+theme_classic()+
  geom_point(data = fieldatabyLoc.gam,aes(fill=market_binomial_data),size=5,
             shape=21,color="black")+sc+xlab("Relative cleaner abundance")+
  ylab("Probability of visitor leaving")+
  labs(fill="Probability \n of choosing \n a visitor")+
  theme(legend.position ="none",
    axis.text = element_text(size=axisSize),
    axis.title.x = element_text(size=axislabSize),
    axis.title.y = element_text(size=axislabSize))

# Panel future reward scatter
scatter.obs.pred.gam<-ggplot(data=fieldatabyLocSamps.gam,
                             aes(y=probVisi.data,colour=site_year,x=probvisitor.pred))+
  stat_gradientinterval(aes(slab_alpha=F),point_interval = mode_hdi,point_size=2)+
  geom_abline(slope=1)+ylab("Observed")+xlab("Predicted")+
  geom_point(data = fieldatabyLoc.gam,aes(y=market_binomial_data,x=probvisitor.pred),
             shape=15,size=2)+
  ylim(0.45,0.85)+xlim(0.45,0.85)+
  guides(color=guide_legend(title="Location"))+
  scale_color_manual(values = multDiscrPallet[1:12])+
  theme_classic()+
  theme(legend.position ="none",#c(.85, .4),#
        axis.text = element_text(size=axisSize),
        axis.title.x = element_text(size=axislabSize),
        axis.title.y = element_text(size=axislabSize),
        legend.key.size = unit(0.20,'cm'))+
  guides(colour=guide_legend(ncol=1))
# ggplot(data = fieldatabyLoc.gam,aes(y=market_binomial_data,x=probvisitor.pred))+
#   # geom_point(data = fieldatabyLocSamps.gam,aes(x=probvisitor.pred,y=probVisi.data),
#   #            color='grey')+
#   stat_halfeye(data = fieldatabyLocSamps.gam,aes(x=probvisitor.pred,y=as.factor(probVisi.data)),
#              color='grey')+
#   geom_point(size=4)+ylim(0.45,0.9)+xlim(0.45,0.9)+
#   geom_abline(slope=1)+ylab("Observed")+xlab("Predicted")+
#   ggtitle("")+
#   theme_classic()+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_text(size=axisSize),axis.title.y = element_text(size=axislabSize),
#         axis.text = element_text(size=14))+
#   geom_text(x = 0.8, y = 0.85, label = expression(y==x), parse = TRUE,size=3)+
#   geom_text(x = 0.8, y = 0.5, label = deparse(bquote(R^2==.(round(rsqr.gam.McFadden,4)))), 
#             parse = TRUE,size=3)

# Panel negative reward contour
cont.obs.pred.Nrew<- ggplot(data = FIAinterpData.Nrew,aes(x=rel.abund.cleaners,y=prob.Vis.Leav,
                                                        fill=market_binomial_data))+
  geom_raster(interpolate = TRUE) +  
  scale_fill_gradientn(limits=c(0.3,1),colours= myPalette(100))+theme_classic()+
  geom_point(data = fieldatabyLoc.gam,aes(fill=market_binomial_data),size=5,
             shape=21,color="black")+sc+xlab("Relative cleaner abundance")+
  ylab("Probability of visitor leaving")+
  labs(fill="Probability \n of choosing \n a visitor")+
  theme(legend.position = "top",
    axis.text = element_text(size=axisSize),
    axis.title.x = element_text(size=axislabSize),
    axis.title.y = element_text(size=axislabSize))

# Panel negative reward scatter


scatter.obs.pred.Nrew<-ggplot(data=fieldatabyLocSamps.Nrew,
                              aes(y=probVisi.data,colour=site_year,x=probvisitor.pred))+
  stat_gradientinterval(aes(slab_alpha=F),point_interval = mode_hdi,point_size=2)+
  geom_abline(slope=1)+ylab("Observed")+xlab("Predicted")+
  geom_point(data = fieldatabyLoc.Nrew,aes(y=market_binomial_data,x=probvisitor.pred),
             shape=15,size=2)+
  ylim(0.45,0.85)+xlim(0.45,0.85)+
  guides(color=guide_legend(title="Location"))+
  scale_color_manual(values = multDiscrPallet[1:12],
                     name = "Location")+
  theme_classic()+theme(legend.position ="top",#c(.9, .4),
                      axis.text = element_text(size=axisSize),
                      axis.title.x = element_text(size=axislabSize),
                      axis.title.y = element_text(size=axislabSize),
                      legend.text = element_text(size=6),
                      legend.key.size = unit(0.1,'cm'))+
  guides(colour=guide_legend(ncol=3,title="",override.aes = list(size=1)))
  
  
# ggplot(data = fieldatabyLoc.Nrew,aes(y=market_binomial_data,x=probvisitor.pred))+
#   geom_point(data = fieldatabyLocSamps.gam,aes(x=probvisitor.pred,y=probVisi.data),
#              color='grey')+
#   geom_point(size=4)+ylim(0.45,0.9)+xlim(0.45,0.9)+
#   geom_abline(slope=1)+ylab("Observed")+xlab("Predicted")+
#   ggtitle("")+
#   theme_classic()+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_text(size=axisSize),axis.title.y = element_text(size=axislabSize),
#         axis.text = element_text(size=axisSize))+
#   geom_text(x = 0.8, y = 0.85, label = expression(y==x), parse = TRUE,size=3)+
#   geom_text(x = 0.8, y = 0.5, label = deparse(bquote(R^2==.(round(rsqr.Nrew.McFadden,4)))), 
#             parse = TRUE,size=3)

# ggarrange(cont.obs.pred.Nrew,scatter.obs.pred.Nrew,
#           labels=c('a','b'),common.legend=FALSE,legend = "top")
# 
# ggarrange(cont.obs.pred.both,scatter.obs.pred.both,
#           cont.obs.pred.gam,scatter.obs.pred.gam,
#           labels=c('a','b','c','d'))

# ggarrange(cont.obs.pred.Nrew,scatter.obs.pred.Nrew,
#           labels=c('a','b'),common.legend=FALSE,legend = "top")
# 
# ggarrange(cont.obs.pred.both,scatter.obs.pred.both,
#           cont.obs.pred.gam,scatter.obs.pred.gam,
#           labels=c('a','b','c','d'))

