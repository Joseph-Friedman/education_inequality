### Joseph Friedman - josephfriedman@mednet.ucla.edu
### Hunter York - hyowrk@uw.edu
##Compile Education Data, Run Cohort Implementation model
########################################################################################
# Setup
########################################################################################
#Install and Load Required Packages
list.of.packages <- c("data.table","stringr","parallel","ggplot2","gridExtra","lme4","splines","boot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
lapply(list.of.packages, require, character.only = TRUE)
library(lmerTest, lib.loc = "<<<<< filepath redacted >>>>>")
library(merTools, lib.loc = "<<<<< filepath redacted >>>>>")

#Arguments
jpath <- "<<<<< filepath redacted >>>>>"
cores <- 10
start_year <- 1950

#Passed arguments
print(commandArgs(trailingOnly = T))
model_version <- commandArgs(trailingOnly = T)[2]
r.dists <- commandArgs(trailingOnly = T)[3]
h.scen <- commandArgs(trailingOnly = T)[4]
c.cohort.mod <- substr(as.character(commandArgs(trailingOnly = T)[5]),1,2)
c.source.adj.mod <- substr(as.numeric(commandArgs(trailingOnly = T)[6]),1,1)
end_year <- substr(as.numeric(commandArgs(trailingOnly = T)[7]),1,4)

c.data.lock <- "none"

#roots
root <- paste0("<<<<< filepath redacted >>>>>")
code_root <- paste0("<<<<< filepath redacted >>>>>")
hroot <- "<<<<< filepath redacted >>>>>"

#load location set
locs <- fread(paste0(root,".csv"))
<<<<< redacted >>>>>
locs<- locs[,.SD,.SDcols=c("iso3","region_name","region_id")]
########################################################################################
#Compile Data
########################################################################################

if (c.data.lock == "none") {
  ##Create List of Files
  ###date of collapsed data to use
  date <- "2019-02-28"
  
  ########################################################################################
  #Compile Education Data
  ########################################################################################
  data.files <- c(list.files(paste0(hroot,"/data/input_data/"),full.names = T,pattern=".csv"),
                  list.files(paste0(jpath,"/binned"),full.names = T),
                  list.files(paste0(jpath,"/tabulated_data/", date, "/binned"),full.names = T),
                  list.files(paste0(jpath,"/tabulated_data/", date),full.names = T, pattern = ".csv"))
  
  #remove extractions that are from 2018 that were redone
  possibledups <- unlist(strsplit(data.files, "/") %>% lapply(tail, n=1))
  dups <-data.files[duplicated(possibledups, fromLast = T)]
  data.files <- data.files[!data.files %in% dups]
  
  ##Load data Using multi-processing if on the cluster
  data.list <- mclapply(data.files,fread,mc.cores=ifelse(Sys.info()[1]=="Windows", 1, cores))
  
  
  #Collapse list of data.tables into one data.table
  edu.dists <- rbindlist(data.list,fill=T)
  
  ########################################################################################
  #Add true zeroes to single-year bins, calculate bigger bins, add means as psuedo-proportion
  ########################################################################################
  #add zeroes
  edu <- edu.dists[,.(proportion=mean(proportion),proportion_se=mean(proportion_se),mean=mean(mean),mean_se=mean(mean_se)),by=.(age_start,sex,year,iso3,sample_size,total_ss,edu_yrs,type,nid)]
  edu[,total_ss := sum(sample_size), by = .(age_start, sex, year, iso3, type, nid)]
  edu[,sum_props := sum(proportion), by = .(age_start, sex, year, iso3, type, nid)]
  edu[,proportion := proportion / sum_props]
  edu[,ts:=paste0(iso3,"&",year,"&",age_start,"&",sex,"&",type, "&", nid, "&", total_ss)] #create unique identifier
  template <- expand.grid(ts=unique(edu$ts),edu_yrs=seq(0,18))
  edu <- merge(edu,template,by=c("ts","edu_yrs"),all=T)
  edu[,ts:=as.character(ts)]
  #clean up missing vars
  edu[is.na(iso3),iso3:=str_split(ts, fixed("&"), n = 7, simplify = T)[,1]]
  edu[is.na(year),year:=str_split(ts, fixed("&"), n = 7, simplify = T)[,2]]
  edu[is.na(age_start),age_start:=str_split(ts, fixed("&"), n = 7, simplify = T)[,3]]
  edu[is.na(sex),sex:=str_split(ts, fixed("&"), n = 7, simplify = T)[,4]]
  edu[is.na(type),type:=str_split(ts, fixed("&"), n = 7, simplify = T)[,5]]
  edu[is.na(nid),nid:=str_split(ts, fixed("&"), n = 7, simplify = T)[,6]]
  edu[is.na(total_ss),total_ss:=str_split(ts, fixed("&"), n = 7, simplify = T)[,7]]
  
  #add 0s for expanded data
  edu[is.na(sample_size),sample_size:=0]
  edu[is.na(proportion),proportion:=0]
  edu[is.na(proportion_se),proportion_se:=0]
  
  #calculate conditional probabilities, add bigger bins
  edu[, total_ss := round(total_ss,3)]
  edu.bins <- dcast.data.table(data=edu[,c("ts","type","iso3","age_start","sex","year","proportion","edu_yrs", "total_ss", "nid")],formula=ts+type+iso3+age_start+sex+year+total_ss+nid~edu_yrs,value.var = "proportion",fun.aggregate = sum)
  #these are used for calculations below
  edu.bins[,prop_1_5:=`1` + `2` +`3` +`4` +`5`]
  edu.bins[,prop_6_11:=`6` + `7` +`8` +`9` +`10` +`11`]
  edu.bins[,prop_12_18:=`12` + `13` +`14` +`15` +`16` +`17`+`18`]
  #entities used to model larger bins
  edu.bins[,cprop_0:=`0`]
  edu.bins[,cprop_12_18:=prop_12_18 / (1-cprop_0)]
  edu.bins[,cprop_1_5:=prop_1_5 / (1-cprop_0+prop_12_18)]
  #calculate conditional probabilities to model single year bins
  #1-5
  edu.bins[,cprop_1:=`1` / prop_1_5]
  edu.bins[,cprop_2:=`2` / (prop_1_5-`1`)]
  edu.bins[,cprop_3:=`3` / (prop_1_5-(`1`+`2`))]
  edu.bins[,cprop_4:=`4` / (prop_1_5-(`1`+`2`+`3`))]
  #6-11
  edu.bins[,cprop_6:=`6` /prop_6_11]
  edu.bins[,cprop_7:=`7` / (prop_6_11-`6`)]
  edu.bins[,cprop_8:=`8` / (prop_6_11-(`6`+`7`))]
  edu.bins[,cprop_9:=`9` / (prop_6_11-(`6`+`7`+`8`))]
  edu.bins[,cprop_10:=`10` / (prop_6_11-(`6`+`7`+`8`+`9`))]
  #12-18
  edu.bins[,cprop_12:=`12` / prop_12_18]
  edu.bins[,cprop_13:=`13` / (prop_12_18-`12`)]
  edu.bins[,cprop_14:=`14` / (prop_12_18-(`12`+`13`))]
  edu.bins[,cprop_15:=`15` / (prop_12_18-(`12`+`13`+`14`))]
  edu.bins[,cprop_16:=`16` / (prop_12_18-(`12`+`13`+`14`+`15`))]
  edu.bins[,cprop_17:=`17` / (prop_12_18-(`12`+`13`+`14`+`15`+`16`))]
  setnames(edu.bins,old=paste0(0:18),new=paste0("prop_",0:18))
  
  #reshape wide
  #edu.bins <- melt.data.table(edu.bins,id.vars=c("ts","type","iso3","age_start","sex","year"),measure.vars = names(edu.bins)[x=grepl(names(edu.bins),pattern="cprop")],variable.name = "code")
  edu.bins <- melt.data.table(edu.bins,id.vars=c("ts","type","iso3","age_start","sex","year", "total_ss", "nid"),variable.name = "code")
  
  
  #add zeroes
  edu.bins[is.na(value),value:=0]

#CREATE means from porportions
  edu.means <- edu[,.(value=mean(mean,na.rm=T)),by=.(age_start,sex,year,iso3,type)]
  edu.means[,code:="mean"]
  
  #transform so that values lie between 0 and 1
  #5-9 year olds max mean is 3
  edu.means[age_start==5 & value>2.95,value:=2.95]
  edu.means[age_start==5,value:=value/3]
  #10-14 year olds max mean is 8
  edu.means[age_start==10 & value>7.95,value:=7.95]
  edu.means[age_start==10,value:=value/8]
  #15-19 year olds max mean is 13
  edu.means[age_start==15 & value>12.95,value:=12.95]
  edu.means[age_start==15,value:=value/13]
  #20+ limit is 18
  edu.means[age_start>15,value:=value/18]
  
  #append data if modeling dists and means
  if (r.dists==1) {
    edu.all <- rbind(edu.bins,edu.means,fill=T)
  } else {
    edu.all <- copy(edu.means)  
  }
  
  #add sample sizes
  ss <- edu[,.(sample_size=sum(sample_size,na.rm=T)),by=.(age_start,sex,year,iso3,type)]
  edu.all <- merge(edu.all,ss,by=c("age_start","sex","year","iso3","type"))
  

  #clean up data provider variable
  edu.all[,type:=gsub(type,pattern="\\/[A-Z][A-Z][A-Z]",replacement="")]
  edu.all[,type:=gsub(type,pattern="\\/[0-9][0-9][0-9]",replacement="")]
  edu.all[,type:=gsub(type,pattern="[A-Z][A-Z][A-Z]\\/",replacement="")]
  edu.all[,type:=gsub(type,pattern="[0-9][0-9][0-9]\\/",replacement="")]
  edu.all[,type:=toupper(type)]
  old <- c("MULTI-COUNTRY SURVEY STUDY ON HEALTH", "PAHO_SABE","SURVEY_ADULT_SKILLS","ISSP","MICS","DHS","IPUMS","MCSS")
  new <- c("MCSS","PAHO_SABE","SURVEY_ADULT_SKILLS","ISSP","MICS","DHS","IPUMS","MCSS")
  for (i in 1:length(old)) {edu.all[grepl(x=type,pattern=old[i]),type:=new[i]]}
  
  
  ########################################################################################
  #Holdouts - first knockout all gold standard for an iso3
  ########################################################################################
  locs <- fread(paste0(root,"ref/locs.csv"))
  <<<<< redacted >>>>>
  locs<- locs[,.SD,.SDcols=c("iso3","region_name","region_id")]
  #1. Gold Standard Knockouts - split DHS and IPUMS into 3 partitions, knockout 1/3 of country-years
  parts <- edu.all[nchar(iso3)==3 & type %in% c("IPUMS","DHS")] #only knockout national data from DHS and IPUMS 
  parts <- merge(parts,locs,by='iso3')
  parts <- unique(parts[,.(iso3, year, type,region_id)])
  parts <- parts[,num:=.N,by=.(type,iso3,region_id)]
  parts <- parts[,.(.N), by = .(iso3, region_id)]
  #stratify over regions to ensure we don't drop all the data
  for (c.reg in unique(parts$region_id)) {parts[region_id==c.reg,part1:=floor((runif(nrow(parts[region_id==c.reg]))*3)+.9999999999)]}
  parts[,N:=NULL]
  parts[,region_id:=NULL]
  edu.all <- merge(edu.all,parts,by=c("iso3"),all=T)
  edu.all[is.na(part1),part1:=0]
  edu.all[!(type %like% "IPUMS") & !(type %like% "DHS"), part1 := 0]
  
  #2. forecasting knockouts - all data post 2011 for locations with data before and after partition
  min.years <- edu.all[,.(min_year=min(year)),by=.(iso3)]
  edu.all <- merge(edu.all,min.years,by='iso3')
  nrow(edu.all[year>min_year & year >= 2012]) / nrow(edu.all) #what proportion of data are we holding out?
  edu.all[year>min_year & year >= 2012,part2:=1]
  edu.all[is.na(part2),part2:=0]
  
  #3. country knockouts - 10% of national locations
  parts <- edu.all[nchar(iso3)==3]
  parts <- merge(parts,locs,by='iso3')
  parts <- parts[,.(num=sum(sample_size)),by=.(iso3,region_id)]
  #stratify over regions to ensure we don't drop all the data
  for (c.reg in unique(parts$region_id)) {parts[region_id==c.reg,part3:=floor((runif(nrow(parts[region_id==c.reg]))*10)+.9999999999)]}
  parts[,num:=NULL]
  parts[,region_id:=NULL]
  edu.all <- merge(edu.all,parts,by=c("iso3"),all=T)
  edu.all[is.na(part3),part3:=0]
  
  #4 cohort testing knockouts - nockout 10% of data with repeat measures across cohort-time
  edu.all[, cohort := year - age_start]
  parts <- edu.all[nchar(iso3)==3]
  parts <- merge(parts,locs,by='iso3')
  parts[age_start > 20, cohort := year - age_start]
  cohort <- parts[!is.na(cohort) & code == "mean",.(num_rep = .N), by = .(iso3, type, cohort, sex)]
  parts <- merge(parts, cohort, by = c("iso3", "type", "cohort", "sex"), all.x = T)
  parts <- parts[,.(num=sum(sample_size)),by=.(type,iso3,year,region_id, sex, num_rep)]
  parts <- parts[num_rep > 1 & type %in% c("DHS", "IPUMS") ]
  parts <- parts[,.(num_rep = max(num_rep)), by = .(type, iso3, year, region_id, sex)]
  #stratify over regions to ensure we don't drop all the data
  for (c.reg in unique(parts$region_id)) {parts[region_id==c.reg,part4:=floor((runif(nrow(parts[region_id==c.reg]))*3 + 9)+.9999999999)]}
  parts[,region_id:=NULL]
  parts[,num_rep := NULL]
  edu.all <- merge(edu.all,parts,by=c("iso3", "type", "year", "sex"),all.x=T)
  edu.all[is.na(part4),part4:=0]
  
  #5 Gold Standard Knockouts - knockout only some of the data from locations with more than one data point
  parts <- edu.all[nchar(iso3)==3 & type %in% c("IPUMS","DHS")] #only knockout national data from DHS and IPUMS 
  parts <- merge(parts,locs,by='iso3')
  parts <- unique(parts[,.(iso3, year, type,region_id)])
  parts <- parts[,num:=.N,by=.(type,iso3,region_id)]
  parts <- parts[,N := .N, by = .(iso3, region_id,type)]
  parts <- parts[N > 1.1]
  parts <- parts[, surv_id := paste0(type, iso3)]
  #stratify over regions to ensure we don't drop all the data
  for (c.reg in unique(parts$surv_id)) {parts[surv_id==c.reg,part5:=floor((runif(nrow(parts[surv_id==c.reg]))*2 + 12)+.9999999999)]}
  #check to ensure not all data from a country got deleted
  redo <- parts[,unique(part5), by = .(iso3,type)] %>% .[,.N, by = .(iso3,type)] %>% .[N == 1, paste0(iso3,type)]
  parts[,maxyear := max(year), by= .(iso3, type)]
  parts[paste0(iso3,type) %in% redo & year < maxyear, part5 := part5 + 1] %>% .[part5 == 15, part5 := 13]
  parts[,N:=NULL]
  parts[,region_id:=NULL]
  parts[,surv_id := NULL]
  parts[,maxyear := NULL]
  parts[,num := NULL]
  edu.all <- merge(edu.all,parts,by=c("iso3", "year", "type"),all=T)
  edu.all[is.na(part5),part5:=0]
  
  #saveRDS(edu.all,paste0(root,"ref/data_lock/",Sys.Date(),".rds")) #save data version set for later use
  
}




########################################################################################
#Remove knockouts
########################################################################################
print(h.scen)
if (h.scen %in% c(1,2,3)) {
  edu.all <- edu.all[part1!=h.scen]  
}  
if (h.scen %in% c(4)) {
  edu.all <- edu.all[part2!=1]  
} 
if (h.scen %in% c(5,6,7,8,9)) {
  edu.all <- edu.all[part3!=h.scen]  
} 
if (h.scen %in% c(10,11,12)) {
  edu.all <- edu.all[part4!=h.scen]  
} 
if (h.scen %in% c(13,14)) {
  edu.all <- edu.all[part5!=h.scen]  
} 

########################################################################################
#Save input data, record number of sources
########################################################################################
saveRDS(edu.all,paste0(hroot,"data/output_data/",model_version,"/input_data/",h.scen,"/raw_input_data.rds"))

#only keep indicators of interest
edu.all <- edu.all[code %in% c("mean", "cprop_0", "prop_0")]


#tabulate country-survey-years
ss <- edu.all[nchar(iso3)==3,.(mean=mean(year,na.rm=T)),by=.(iso3, year,type)]
ss[,obs:=  1]
ss <- ss[,.(Number=sum(obs)),by=.(type)]
#write.csv(ss[order(-Number)],paste0(hroot,"data/output_data/",model_version,"/input_data/",h.scen,"/sources_record.csv") )

########################################################################################
#Only Keep Modeled Quantities 
########################################################################################

# merge on numeric codes for quantities 0-18
edu.codes <- fread(paste0(root,"ref/edu_codes.csv"))
edu.all <- merge(edu.all,edu.codes,by=c("code"))
edu.all <- edu.all[,.SD,.SDcols=c("iso3","year","age_start","sex","type","edu_yrs","sample_size","value","part1","part2","part3", "part4", "part5")]
edu.all <- edu.all[order(iso3,year,age_start,sex,type,edu_yrs)]


########################################################################################
#Cohort Imputation -- Draws on module files, creates edu.all.imp object
########################################################################################
source(paste0(code_root, "modules/cohort_change2.R"), echo = T) #model cohort change


#######################################################################################
#Source Adjustment Model
#######################################################################################
if(c.source.adj.mod ==4){
  source(paste0(code_root, "modules/type_adjust_model.R"),echo=T) #explore differences by cohort
}
########################################################################################
#merge on location information, save prepped data
########################################################################################
# locs <- get_location_metadata(location_set_id = 22)
locs <- fread(paste0(root,"ref/locs.csv"))
<<<<< redacted >>>>>
locs<- locs[,.SD,.SDcols=c("iso3","region_name","region_id")]
edu.all.imp <- data.table(merge(edu.all.imp,locs))
edu.all.imp[,reg_sex:=paste0(region_id,"_",sex)]
edu.all.imp[,iso3_age:=paste0(iso3,age_start)]



#save region-sex specific files
savefiles <-function(c.reg_sex){saveRDS(edu.all.imp[reg_sex==c.reg_sex],paste0(hroot,"data/output_data/",model_version,"/input_data/",h.scen,"/",c.reg_sex,".rds") )}
mclapply(unique(edu.all.imp$reg_sex), savefiles, mc.cores = cores)