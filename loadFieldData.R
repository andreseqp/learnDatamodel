## Merge field and experimental data to use in the bayesian analysis

library("data.table")
library("here")
library("readxl")
library('dplyr')
library('ggplot2')

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
fieldData.cleaner<-data.table(read_xlsx(here("Data","market_raw_data.xlsx"),sheet = "round_data"))

fieldData.cleaner[,unique(site_year)]

# Fix mistake in location code name
fieldData.cleaner[site_year=="NHS2017",site_year:="NHS 2017"]



# Choose data from the last 20 experimental choices
fieldData.cleaner[,use_sims:=ifelse(is.na(use_sims),0,1)]
fieldData.cleaner[,use_sims:=as.logical(use_sims)]

fieldData.filt<-fieldData.cleaner[use_sims==T,]
fieldData.filt[,score_visitor:=as.numeric(score_visitor)]
fieldData.filt[,stage:=as.factor(stage)]


str(fieldData.filt)

# Sum choices for each location and experimental subject
fieldData.sum<-fieldData.filt[,sum(score_visitor),
                              by=.(site_year,cleaner_ID)]

str(fieldData.sum)

# Rename columns for convenience
names(fieldData.sum)[3]<-"score_visitor"
names(fieldData.sum)[1]<-"site.year"

fieldData.sum[,unique(site.year)]
fieldData.site$site_year



# Get environmental variables from field data
fieldData.sum[,abund.cleaners:=
                fieldData.site[match(site.year,site_year),abundance_cleaners_100m2]]

# Check the location names coincide
fieldData.sum[order(abund.cleaners),max(abund.cleaners),by=site.year]$site.year==
  
  fieldData.site[order(abundance_cleaners_100m2),.(site_year,abundance_cleaners_100m2)]$site_year

fieldData.sum[,abund.visitors:=
                fieldData.site[match(site.year,site_year),abundance_large_100m2]]

fieldData.sum[,abund.residents:=
                fieldData.site[match(site.year,site_year),abundance_small_100m2]]

fieldData.sum[,prob.Vis.Leav:=
                fieldData.site[match(site.year,site_year),percentage_swim_off*0.01]]
# Standardize location and cleaner IDs
fieldData.sum[,site.year:=gsub(" ",replacement = "_",x = site.year)]
fieldData.sum[,cleaner_ID:=gsub(" ",replacement = "",x = cleaner_ID)]

str(fieldData.sum)

# Organize cleaners on competent and incompetent
threashold<-1.5
fieldData.sum[,highDens:=abund.cleaners>threashold]
fieldData.sum[,visPref:=sapply(score_visitor,function(x){
  binom.test(x,n=20,p=0.5,alternative = "greater")$p.value<0.05
})]
fieldData.sum[,competence:=highDens+visPref!=1]


setorder(fieldData.sum,site.year)

fieldData.sum[,`:=`(highDens=NULL,visPref=NULL,competence=as.integer(competence))]

# Print dataset with both sources of information with competence grouping
fwrite(fieldData.sum,here("data",paste0("data_cleaner_abs_threa",threashold,".txt")),
       row.names = FALSE,sep = "\t")

# Print dataset with both sources of information without competence grouping
fieldData.sum[,competence:=NULL]
fwrite(fieldData.sum,here("data","data_ABC_cleaner_absolute.txt"),
       row.names = FALSE,sep = "\t")

## 
