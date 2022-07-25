## Merge field and experimental data to use in the bayesian analysis

library("data.table")
library("here")
library("readxl")
library('dplyr')
library('ggplot2')
library('cowplot')
library("ggbeeswarm")
source(here("aesth_par.R"))

## Load field data based on the market experiment performance criteria
fieldData<-fread(here("Data","market_model.csv"))


str(fieldData)
## Summarize ecological data by site

fieldData[,unique(site_year)]
names(fieldData)

# Summarize field data by location
fieldData.site<-fieldData[,lapply(.SD,max),by=site_year,
                          .SDcols=c("abundance_large_100m2","abundance_small_100m2",
                                   "abundance_cleaners_100m2","percentage_swim_off")]


## Load the experimental data
fieldData.cleaner<-data.table(read_xlsx(here("Data","market_raw_data.xlsx"),
                                        sheet = "round_data"))

fieldData.cleaner[,unique(site_year)]

# Fix mistake in location code name
fieldData.cleaner[site_year=="NHS2017",site_year:="NHS 2017"]

# standarize cleaner IDs and reef-year sites
fieldData.cleaner[,cleaner_ID:=gsub(" ",replacement = "",x = cleaner_ID)]
fieldData.cleaner[,cleaner_ID:=gsub("c",replacement = "C",x = cleaner_ID)]
fieldData.cleaner[,cleaner_ID:=gsub("l",replacement = "L",x = cleaner_ID)]
fieldData.cleaner[,site_year:=gsub(" ",replacement = "_",x = site_year)]

# get scores from last two session, both initial and reversal
fieldData.cleaner.filt<-fieldData.cleaner[,.SD[c(.N-1,.N)],
                                          by=.(site_year,cleaner_ID,
                                               stage)]
fieldData.cleaner.filt<-fieldData.cleaner.filt[stage!="NA"]
fieldData.cleaner.filt[,score_visitor:=as.numeric(score_visitor)]


fieldData.cleaner.filt.sum<-fieldData.cleaner.filt[,length(score_visitor),
                       by=.(site_year,cleaner_ID)]

fieldData.cleaner.filt.sum[V1>2,length(V1)]+
  fieldData.cleaner.filt.sum[V1==2,length(V1)]
# get the two sessions for analysis based on rules
# last two sessions if initial was not solved
# Last from initial and reversal if initial solved
# get length of the vector to check
# classify performance: 
#       0 initial not solved
#       1 initial solved but not reversal
#       2 both initial and reversal solved
fieldData.cleaner.sum<-fieldData.cleaner.filt[,
.(score_visitor=
    ifelse(sum(stage=="reversal")<1,sum(score_visitor),0)+
    ifelse(sum(stage=="reversal")>1&length(stage)==4,
           sum(score_visitor[c(2,4)]),0)+
    ifelse(sum(stage=="reversal")==1&length(stage)==3,
           sum(score_visitor[c(2,3)]),0)+
    ifelse(sum(stage=="reversal")==2&length(stage)==3,
           sum(score_visitor[c(1,3)]),0),
  class.perform=
    ifelse(sum(stage=="reversal")<1,0,
      ifelse(sum(stage=="reversal")>=1&max(round)<20,2,1)),
             l.score=length(score_visitor))
                                          ,by=.(site_year,cleaner_ID)]

fieldData.cleaner.sum

# Choose data from the last session manually
fieldData.cleaner[,use_sims:=ifelse(use_sims=="NA",FALSE,TRUE)]

fieldData.filt.man<-fieldData.cleaner[use_sims==T,]
fieldData.filt.man[,score_visitor:=as.numeric(score_visitor)]
fieldData.filt.man[,stage:=as.factor(stage)]


str(fieldData.filt.man)

# Sum choices for each location and experimental subject
fieldData.sum.man<-fieldData.filt.man[,.(score_visitor.man=sum(score_visitor)),
                              by=.(site_year,cleaner_ID)]

str(fieldData.sum.man)

fieldData.cleaner.sum<-fieldData.cleaner.sum[fieldData.sum.man[,.(site_year,cleaner_ID,score_visitor.man)],
                                                       on=.(site_year=site_year,cleaner_ID=cleaner_ID)]

# check that manual and code coincide
fieldData.cleaner.sum[score_visitor!=score_visitor.man]


# Rename columns for convenience
names(fieldData.cleaner.sum)[1]<-"site.year"

fieldData.cleaner.sum[,unique(site.year)]

fieldData.site[,site_year:=gsub(" ",replacement = "_",x = site_year)]

# Get environmental variables from field data
fieldData.cleaner.sum[,abund.cleaners:=
                fieldData.site[match(site.year,site_year),abundance_cleaners_100m2]]

# Check the location names coincide
fieldData.cleaner.sum[order(abund.cleaners),max(abund.cleaners),by=site.year]$site.year==
  
  fieldData.site[order(abundance_cleaners_100m2),.(site_year,abundance_cleaners_100m2)]$site_year

fieldData.cleaner.sum[,abund.visitors:=
                fieldData.site[match(site.year,site_year),abundance_large_100m2]]

fieldData.cleaner.sum[,abund.residents:=
                fieldData.site[match(site.year,site_year),abundance_small_100m2]]

fieldData.cleaner.sum[,prob.Vis.Leav:=
                fieldData.site[match(site.year,site_year),percentage_swim_off*0.01]]


str(fieldData.cleaner.sum)

# Organize cleaners on competent and incompetent
threashold<-1.5
fieldData.cleaner.sum[,highDens:=abund.cleaners>threashold]
fieldData.cleaner.sum[,visPref:=sapply(score_visitor,function(x){
  binom.test(x,n=20,p=0.5,alternative = "greater")$p.value<0.05
})]
fieldData.cleaner.sum[,competence:=highDens+visPref!=1]


setorder(fieldData.cleaner.sum,site.year)

# fieldData.cleaner.sum.print<-fieldData.cleaner.sum[,
#           `:=`(highDens=NULL,visPref=NULL,
#                class.perform=NULL,l.score=NULL,
#                score_visitor.man=NULL,
#                competence=as.integer(competence))]


# Print dataset with both sources of information with competence grouping
# fwrite(fieldData.cleaner.sum.print,here("data",paste0("data_cleaner_abs_threa",threashold,".txt")),
       # row.names = FALSE,sep = "\t")

# Print dataset with both sources of information without competence grouping
# fieldData.cleaner.sum.print[,competence:=NULL]
# fwrite(fieldData.cleaner.sum.print,here("data","data_cleaner_absolute.txt"),
       # row.names = FALSE,sep = "\t")

## Recalculate preference for visitors based only on the initial round ---------

fieldData.initial<-fieldData.cleaner[stage=="initial"]

fieldData.initial.last2<-fieldData.initial[,.SD[c(.N-1,.N)],
                                           by=.(site_year,cleaner_ID)]

countFieldata<-fieldData.initial.last2[,length(score_visitor),by=.(site_year,cleaner_ID)]

fieldData.initial.last2[,score_visitor:=as.numeric(score_visitor)]

fieldData.in.last2.sum<-fieldData.initial.last2[,.(score_visitor.init=sum(score_visitor),
                                                   nRounds.init=length(score_visitor)),
                                                by=.(site_year,cleaner_ID)]

fieldData.cleaner.sum[,unique(site.year)]
fieldData.in.last2.sum[,unique(site_year)]

fieldData.merged<-fieldData.cleaner.sum[fieldData.in.last2.sum,on=.(site.year=site_year,
                                          cleaner_ID=cleaner_ID)]

fieldData[,site_year:=gsub(" ",replacement = "_",x = site_year)]
fieldData[,ID:=gsub(" ",replacement = "",x = ID)]
fieldData[,ID:=gsub("c",replacement = "C",x = ID)]
fieldData[,ID:=gsub("l",replacement = "L",x = ID)]



field.merged.2<-fieldData.merged[fieldData[,.(site_year,ID,market_binomial,market_both_score)],
                 on=.(site.year=site_year,cleaner_ID=ID)]

# png(here("FigsExtra","initVSrev.png"),width = 500,height = 500)

ggplot(field.merged.2,aes(x=score_visitor.init,y=score_visitor,
                            color=as.factor(market_binomial),
                            shape=as.factor(nRounds.init)))+
  geom_point(size=2)+theme_classic()+
  xlab("score using only initial")+ylab("Score using initial and reversal")+
  labs(color="Performance in \n the experiment",
       shape="number of initial \n rounds available")+
  scale_colour_manual(values = multDiscrPallet[1:2])

# dev.off()


# png(here("FigsExtra","scoresVSperform.png"),width = 800,height = 800)

initANDRev<-ggplot(field.merged.2,aes(y=score_visitor,x=as.factor(class.perform)))+
  geom_beeswarm()+theme_classic()+labs(title = "Initial and reversal")+
  scale_x_discrete(name ="Performance",labels=c("0"="None","1"="initial only",
                                                  "2"="both"))+ylim(0,20)+
  stat_summary(fun.data = function(x){return(c(y=20,
               label=length(x)))},
               geom = "text")+ylab("# of visitor chosen")+
  theme(plot.title = element_text(hjust = 1)) 


init<-ggplot(field.merged.2,aes(y=score_visitor.init,x=as.factor(class.perform)))+
  geom_beeswarm()+theme_classic()+labs(title = "Initial")+
  scale_x_discrete(name ="Performance",labels=c("0"="None",
                                                "1"="initial only",
                                                "2"="both"))+ylim(0,20)+
  stat_summary(fun.data = function(x){return(c(y=20,
                          label=length(x)))},
               geom = "text")+ylab("")+
  theme(plot.title = element_text(hjust = 1)) 



# plot_grid(initANDRev,init,ncol = 2)
# dev.off()
