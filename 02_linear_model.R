### Joseph Friedman - josephfriedman@mednet.ucla.edu
### Hunter York - hyowrk@uw.edu
##Compile Education Data, Fit linear model 
#Install and Load Required Packages
list.of.packages <- c("data.table","stringr","gridExtra","boot","lme4","splines","ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)
#Arguments
jpath <- "<<<<< filepath redacted >>>>>"
cores <- 20
start_year <- 1950
root <- paste0("<<<<< filepath redacted >>>>>")
code_root <- paste0(root, "code")
hroot <- "<<<<< filepath redacted >>>>>"



#Passed arguments
print(commandArgs())
model_version <- commandArgs()[3]
c.reg_sex <- commandArgs()[4]
h.scen <- commandArgs()[5]
end_year <- as.numeric(commandArgs()[6])

#load data
edu.orig <- readRDS(paste0(hroot,"data/output_data/",model_version,"/input_data/",h.scen,"/",c.reg_sex,".rds"))

#holdout scenarios (repetitive with part 1, but just in case)
print(h.scen)
if (h.scen %in% c(1,2,3)) {
  edu.dat <- edu.orig[part1!=h.scen]  
}   
if (h.scen %in% c(4)) {
  edu.dat <- edu.orig[part2!=1]  
} 
if (h.scen %in% c(5,6,7,8,9)) {
  edu.dat <- edu.orig[part3!=h.scen]  
} 
if (h.scen %in% c(0)) {
  edu.dat <- copy(edu.orig)
} 
if (h.scen %in% c(10,11,12)) {
  edu.dat <- edu.orig[part4!=h.scen]  
} 
if (h.scen %in% c(13,14)) {
  edu.dat <- edu.orig[part5!=h.scen]  
} 
##########################################################################################################
#Merge onto prediction template
##########################################################################################################
#locs <- get_location_metadata(location_set_id = 22)
locs <- fread(paste0(root,"ref/locs.csv"))
locs <- locs[level>2]
<<<<< redacted for peer review >>>>>
locs<- locs[,.SD,.SDcols=c("iso3","region_name","region_id","level")]
locs <- locs[region_id==as.numeric(strsplit(c.reg_sex,split="_")[[1]][1])]

template <- expand.grid(iso3=unique(locs$iso3),year=seq(start_year,end_year,1),age_start=unique(edu.dat$age_start),sex=unique(edu.dat$sex),edu_yrs=unique(edu.dat$edu_yrs))
edu.dat <- merge(edu.dat,template,by=c("age_start","sex","year","iso3","edu_yrs"),all.y=T)

edu.dat[,reg_sex:=paste0(region_id,"_",sex)]
edu.dat[,iso3_age:=paste0(iso3,age_start)]
edu.dat <- edu.dat[order(iso3,year,age_start,sex)]

##########################################################################################################
#Model - edu_year == 18 is modeling means, otherwise modeling proportions
##########################################################################################################
#proportions
edu.dat[value>.999,value:=.999]
edu.dat[value<.001,value:=.001]
edu.dat[,logit_value:=logit(value)]
edu.dat[,iso3:=as.character(iso3)]


options(lmerControl=list(check.nobs.vs.rankZ = "warning",
                         check.nobs.vs.nlev = "warning",
                         check.nobs.vs.nRE = "warning",
                         check.nlev.gtreq.5 = "warning",
                         check.nlev.gtr.1 = "warning"))

for (c.edu_year in unique(edu.dat$edu_yrs)) {
  print(c.edu_year)
  if (c.edu_year== 18) {
    mod <- lmer(logit_value~ year + (1|iso3) +ns(age_start,df=edu.dat[edu_yrs == c.edu_year],knots=c(45,65)),data=edu.dat[edu_yrs == c.edu_year])
    edu.dat[edu_yrs == c.edu_year,logit_value_pred:= predict(mod,newdata=edu.dat[edu_yrs==c.edu_year],allow.new.levels=T)]
  } else {
    mod <- try(lmer(logit_value~ year + (1|iso3_age) +ns(age_start,df=edu.dat[edu_yrs==c.edu_year],knots=c(45,65)),data=edu.dat[edu_yrs==c.edu_year]))
    if (class(mod) == "try-error"){mod <- lmer(logit_value~ year + (1|iso3_age),data=edu.dat[edu_yrs==c.edu_year])}
    edu.dat[edu_yrs == c.edu_year,logit_value_pred:= predict(mod,newdata=edu.dat[edu_yrs==c.edu_year],allow.new.levels=T)]
  }
} 

#calculate region-sex MAD for GPR
reg_mad <- copy(edu.dat)
reg_mad[,abs.dev := abs(logit_value - logit_value_pred)]
reg_mad = reg_mad[,.(mad=median(abs.dev,na.rm=T)),by=(edu_yrs)]
#edu[,mad:=reg_mad$mad[1]]
edu.dat <- merge(edu.dat,reg_mad,by='edu_yrs')
edu.dat[,mad:=mad]

#calculate and delta transform variance
#proportions
edu.dat[,variance:=(((value*(1-value))/sample_size)^2)]
edu.dat[,extrap_var := extrap_ci^2]
edu.dat[,variance := variance + extrap_var]


##add in extrapolation variance

edu.dat[,logit_variance:=variance * (1/((value)*(1-(value))))^2+mad]
edu.dat[value < .01,logit_variance:=variance * (1/((.01)*(1-(.01))))^2+mad]
#means
edu.dat[,mn_variance:=((value*(1-value))/sample_size)^2]
edu.dat[,mn_variance := mn_variance + extrap_var]

edu.dat[,mn_logit_variance:=mn_variance * (1/((value)*(1-(value))))^2+mad]
edu.dat[value < .01,mn_logit_variance:=mn_variance * (1/((.01)*(1-(.01))))^2+mad]
edu.dat[edu_yrs==18,logit_variance:=mn_logit_variance]
edu.dat[edu_yrs==18,variance:=mn_variance]

#remove vars
edu.dat[,extrap_var := NULL]
edu.dat[,extrap_ci := NULL]

###standardize panel variables
<<<<< redacted for peer review >>>>>
  edu.dat[age_start < 80,age_group_id:=((age_start - 15) /5 ) + 8]
edu.dat[age_start > 75,age_group_id:=((age_start - 15) /5 ) + 17]
edu.dat[age_start ==95,age_group_id:=235]

# locs <- get_location_metadata(location_set_id = 22)
locs <- fread(paste0(root,"ref/locs.csv"))
<<<<< redacted for peer review >>>>>
  <<<<< redacted for peer review >>>>>
  
save.data[,reg_sex:=NULL]
save.data[,region_id:=NULL]
save.data[,region_name:=NULL]
save.data[,iso3_age:=NULL]
save.data[,mn_variance:=NULL]
save.data[,mn_logit_variance:=NULL]
save.data <- save.data[order(<<<<< redacted for peer review >>>>>,age_start,sex_id,edu_yrs,year_id)]


var05th <- quantile(save.data$logit_variance,probs=.05,na.rm=T)
save.data[logit_variance < var05th,logit_variance:= var05th ]
var95th <- quantile(save.data$logit_variance,probs=.95,na.rm=T)
save.data[logit_variance > var95th,logit_variance:= var95th ]

# some dups exist bc saving data too not just preds
save.data <- save.data[age_group_id %in%  c(6:20, 30:32, 235)]
#assert_ids(save.data, id_vars = list(<<<<< redacted for peer review >>>>> = template$iso3, age_group_id = c(6:20, 30:32, 235), sex_id = 1:2, year_id = 1950:2040), assert_dups = F)

save.data$sex_id <- as.integer(save.data$sex_id)
save.data$sample_size <- as.integer(save.data$sample_size)
save.data$sex_id <- as.integer(save.data$sex_id)


#Save region-sex specific preds; launch GPR
dir.create(paste0(hroot,"data/output_data/",model_version,"/linear/",h.scen,"/"),showWarnings = F,recursive = T)
dir.create(paste0(hroot,"data/output_data/",model_version,"/gpr/",h.scen,"/"),showWarnings = F,recursive = T)

write.csv(save.data,paste0(hroot,"data/output_data/",model_version,"/linear/",h.scen,"/",c.reg_sex,".csv"),row.names=F)




