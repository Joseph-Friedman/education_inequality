## Process - create single-year distributions from modelled means and proportion with 0 years of education
### Joseph Friedman - josephfriedman@mednet.ucla.edu
### Hunter York - hyowrk@uw.edu
#Install and Load Required Packages
list.of.packages <- c("data.table","stringr","ggplot2","gridExtra","boot", "parallel", "reldist", "assertable")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
lapply(list.of.packages, require, character.only = TRUE)
#Arguments
jpath <- "<<<<< filepath redacted >>>>>"

library(merTools, lib.loc = "/<<<<< filepath redacted >>>>>/")

start_year <- 1950
#Passed arguments
model_version <- commandArgs(trailingOnly = T)[2]
c.draws <- as.numeric(commandArgs(trailingOnly = T)[3])
h.scen <- as.numeric(commandArgs(trailingOnly = T)[4])
r.rake <- as.numeric(commandArgs(trailingOnly = T)[5])
forecasting <- substr(commandArgs(trailingOnly = T)[6],1,4)
end_year <-  as.numeric(substr(commandArgs(trailingOnly = T)[7],1,4))

task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
parameters <- fread("/<<<<< filepath redacted >>>>>/gpr_parameters.csv")

c.loc_sex <- parameters[task_id, loc_sex]

#create holdouts for forecasting if necessary
if(h.scen == 4){year_before_forecasts <- 2008}else{year_before_forecasts <- 2019}

root <- paste0(<<<<< filepath redacted >>>>>)
code_root <- paste0("<<<<< filepath redacted >>>>>", "/code")
hroot <- "/<<<<< filepath redacted >>>>>/"

#############################################################################################################################################################
#load data
#############################################################################################################################################################
print(paste0(hroot,"/data/output_data/",model_version,"/gpr/",h.scen,"/",c.loc_sex,".csv"))
print(c(model_version, c.draws, h.scen, r.rake, forecasting))

#load data from GPR steps
  edu <- fread(paste0(hroot,"/data/output_data/",model_version,"/gpr/",h.scen,"/",c.loc_sex,".csv"))
  
  edu$V1 <-NULL
  
  #reshape long if using draws
  if("draw" %in% substr(names(edu),1,4)){
    edu <- melt(edu, id.vars = c("<<<redacted>>>", "edu_yrs", "age_start", "sex_id", "year_id", "type", "logit_value_pred"), measure.vars = patterns("draw"))
    edu[, draw := variable] %>% .[,variable := NULL] %>% .[,gpr_mean := value] %>% .[,value := NULL]
    cores <- 10
  }else {
    edu$draw <- "draw0"
    cores <- 5
  }
  
  edu.dat <- edu[,.(gpr_mean=mean(gpr_mean),linear=mean(logit_value_pred)),by=.(<<<redacted>>>,year_id,age_start,sex_id,edu_yrs,draw)]
  #merge on meaning
  edu.codes <- fread(paste0(root,"ref/edu_codes.csv"))
  edu.dat <- merge(edu.dat,edu.codes,by=c("edu_yrs"), all.x = T)
  
  edu.dat <- melt.data.table(edu.dat,id.vars=c("<<<redacted>>>","year_id","age_start","sex_id","edu_yrs","code","description", "draw"),variable.name = "method")
  
  #get the out of logit space
  edu.dat[,value:=inv.logit(value)]
  
  #Create proportions from conditional probabilities 
  edu.dat <- dcast.data.table(data=edu.dat,formula=<<<redacted>>>+age_start+sex_id+year_id+method+draw~code,value.var = "value")

#############################################################################################################################################################
#process to regular proportions and means
#############################################################################################################################################################

#enforce monotonocity 0 years of education decreasing
edu.dat <- edu.dat[order(method,<<<redacted>>>,age_start,sex_id,year_id,draw)]

#larger bins
edu.dat[,prop_0:=cprop_0]
edu.dat <- edu.dat[method == "gpr_mean"]
mins<- edu.dat[,minimum := min(prop_0), by = .(age_start, sex_id, <<<redacted>>>, draw, method)]
min_years <- edu.dat[prop_0 == minimum, .(year_id, prop_0), by = .(age_start, sex_id, <<<redacted>>>, draw, method)]
min_years <- min_years[,.(prop_0 = min(prop_0), year_id = min(year_id)), by = .(age_start, sex_id, <<<redacted>>>, draw, method)]
setnames(min_years, c("prop_0", "year_id"), c("new_prop_0", "min_year"))
edu.dat <- merge(min_years, edu.dat, by = c("age_start", "sex_id", "<<<redacted>>>", "draw", "method"), all.y = T)
edu.dat[year_id > min_year, prop_0 := new_prop_0]
edu.dat[, c("minimum", "new_prop_0", "min_year") := NULL]
edu.dat[, cprop_0 := prop_0]


##launch knn algorithm
source(paste0(code_root, "modules/knn_prop0.R" ), echo = T)



#transform means to correct space using ceilings 
#5-9 year olds max mean is 3
edu.dat[age_start==5,mean:=mean*3]
#10-14 year olds max mean is 8
edu.dat[age_start==10,mean:=mean*8]
#15-19 year olds max mean is 13
edu.dat[age_start==15,mean:=mean*13]
#20+ limit is 18
edu.dat[age_start>15,mean:=mean*18]

#melt long for saving
edu.dat[,l1d:=NULL]
edu.dat.wide <- edu.dat[,.SD, .SDcols = c("year_id", "sex_id", "<<<redacted>>>", "age_start", "method", "draw", names(edu.dat)[names(edu.dat) %like% "prop"], "mean", "weighted_mean")]
edu.data.long <- melt.data.table(edu.dat.wide,id.vars=c("method","<<<redacted>>>","age_start","sex_id","year_id","draw"))
edu.codes[,variable:=code]
edu.data.long <- merge(edu.data.long,edu.codes,by=c("variable"),all.x = T)


#####forecast using ROC. Only start forecasting for cohorts which have not yet completed education (<25 years old in 2019)
edu.data.long <- edu.data.long[(year_id <= year_before_forecasts & age_start <= 25)|
                                 (year_id <= year_before_forecasts + 5 & age_start == 30)|
                                 (year_id <= year_before_forecasts + 10 & age_start == 35)|
                                 (year_id <= year_before_forecasts + 15 & age_start == 40)|
                                 (year_id <= year_before_forecasts + 20& age_start == 45)|
                                 (year_id <= year_before_forecasts + 25& age_start == 50)|
                                 (year_id <= year_before_forecasts + 30& age_start == 55)|
                                 (year_id <= year_before_forecasts + 35& age_start == 60)|
                                 (year_id <= year_before_forecasts + 40& age_start == 65)|
                                 (year_id <= year_before_forecasts + 45& age_start == 70)|
                                 (year_id <= year_before_forecasts + 50& age_start == 75)|
                                 (year_id <= year_before_forecasts + 55& age_start == 80)|
                                 (year_id <= year_before_forecasts + 60& age_start == 85)|
                                 (year_id <= year_before_forecasts + 65& age_start == 90)|
                                 (year_id <= year_before_forecasts + 70& age_start == 95)]
edu.data.long[variable %in% paste0("prop_",c(0:18, "1_5", "6_11", "12_18")) & value < .001, value := .001]
edu.data.long[variable %in% paste0("prop_",c(0:18, "1_5", "6_11", "12_18")), logit_value := as.numeric(logit(value))]

#create lagged variable to create instantanous rate of change
edu.data.long[, l1d := data.table::shift(logit_value, 1), by = .(variable, method, <<<redacted>>>, age_start, sex_id,draw)]
edu.data.long[, roc := logit_value - l1d]

#only start forecasting when we run out of data
year_bounds <- edu.data.long[,.(max = max(year_id)), by = .(age_start)]
edu.data.long <- merge(edu.data.long, year_bounds, by = "age_start")

#use median rate of change for last 15 years to predict into the future
avg.rocs <- edu.data.long[,.(avg.roc = median(roc[year_id %in% (max - 15):(max - 1)], na.rm = T), value_2018 = logit_value[year_id == max]), by=.(draw,variable, method, <<<redacted>>>, age_start, sex_id, max) ]

#create template to predict on
forecasts <- expand.grid(variable = unique(edu.data.long$variable),
                         method = "gpr_mean",
                         <<<redacted>>> = unique(edu.data.long$<<<redacted>>>),
                         age_start = unique(edu.data.long$age_start),
                         year_id = (year_before_forecasts+1):end_year,
                         sex_id = unique(edu.data.long$sex_id),
                         draw = unique(edu.data.long$draw)) %>% 
  as.data.table()

#only keep age-years which are of interest (which we can't infer from cohort extrapolation)
forecasts <- forecasts[(year_id > year_before_forecasts & age_start <= 25)|
                         (year_id > year_before_forecasts + 5 & age_start == 30)|
                         (year_id > year_before_forecasts + 10 & age_start == 35)|
                         (year_id > year_before_forecasts + 15& age_start == 40)|
                         (year_id > year_before_forecasts + 20& age_start == 45)|
                         (year_id > year_before_forecasts + 25& age_start == 50)|
                         (year_id > year_before_forecasts + 30& age_start == 55)|
                         (year_id > year_before_forecasts + 35& age_start == 60)|
                         (year_id > year_before_forecasts + 40& age_start == 65)|
                         (year_id >  year_before_forecasts + 45& age_start == 70)|
                         (year_id >  year_before_forecasts + 50& age_start == 75)|
                         (year_id >  year_before_forecasts + 55& age_start == 80)|
                         (year_id >  year_before_forecasts + 60& age_start == 85)|
                         (year_id >  year_before_forecasts + 65& age_start == 90)|
                         (year_id >  year_before_forecasts + 70& age_start == 95)]

forecasts <- merge(forecasts, avg.rocs, by = c("variable", "method", "<<<redacted>>>", "age_start", "sex_id","draw"), all.x = T)
forecasts[is.nan(avg.roc), avg.roc := 0]

#apply roc using value in 2019 (coded as value_2018)
  forecasts[variable %in% paste0("prop_",c(0:18, "1_5", "6_11", "12_18")), logit_value_pred := value_2018 + ((avg.roc) * (year_id - max))]
  forecasts[, value := inv.logit(logit_value_pred)]
  forecasts[,max := NULL]
  
  #rake twice, first to chunky then to granular (ensure all larger bins add to 1)
  forecasts[variable %in%paste0("prop_",c(0, "1_5", "6_11", "12_18")),  rake_val1 := sum(value[variable %in%paste0("prop_",c("1_5", "6_11", "12_18"))], na.rm = T)/(1-value[variable == "prop_0"]), by= .(method, <<<redacted>>>, age_start, sex_id, year_id, draw)]
  forecasts[variable %in%paste0("prop_",c("1_5", "6_11", "12_18")), value := value /rake_val1]
  forecasts[variable %in% paste0("prop_",c(0)), rake_var := "prop_0"]
  forecasts[variable %in% paste0("prop_",c(1:5)), rake_var := "prop_1_5"]
  forecasts[variable %in% paste0("prop_",c(6:11)), rake_var := "prop_6_11"]
  forecasts[variable %in% paste0("prop_",c(12:18)), rake_var := "prop_12_18"]
  
  forecasts[,sub_tot := sum(value), by = .(method, <<<redacted>>>, age_start, sex_id, year_id, draw, rake_var)]
  targets <- forecasts[variable %in%paste0("prop_",c(0, "1_5", "6_11", "12_18")), .(variable, method, <<<redacted>>>, age_start, sex_id, draw,year_id, value)]
  setnames(targets, c("value", "variable"), c("target", "rake_var"))
  forecasts <- merge(targets, forecasts, all.y = T, by = c("rake_var", "method", "<<<redacted>>>", "age_start", "sex_id", "draw", "year_id"))
  forecasts[variable %in% paste0("prop_",0:18) & is.na(value), value:= .0001]
  forecasts[variable %in% paste0("prop_",0:18), rake_val3 := sum(value[variable %in% paste0("prop_",1:18)])/(1-value[variable == "prop_0"]), by= .(method, <<<redacted>>>, age_start, sex_id, year_id, draw)]
  #rake a second time to ensure smaller bins add up to larger bins
  forecasts[variable %in% paste0("prop_",1:18), value := value /rake_val3]
  forecasts[, rake_val1 := NULL]
  forecasts[, rake_val3 := NULL]
  forecasts[, sub_tot := NULL]
  forecasts[, target := NULL]
  forecasts[, rake_var := NULL]
  forecasts[, logit_value_pred := NULL]
  forecasts[, logit_value := NULL]


#exclude forecasted summary variables and calculate these from proportions (ergo dont forecast mean directly, instead construct it from forecasted bins)
forecasts <- forecasts[!variable%in% c("mean", "weighted_mean")]
forecasts[,avg.roc := NULL][,value_2018:= NULL]
edu.data.long[, l1d := NULL][,roc:= NULL]
edu.data.long<- rbind(edu.data.long, forecasts, fill = T)
edu.data.long[is.na(code), code := variable]


#set logit floors for small values
edu.data.long[value <= .003, value := .0001]

#copy values to variables
edu.data.long[variable == "mean", value := mean]

#a final sanity check to ensure raking worked
edu.data.long[variable %in% paste0("prop_",0:18), rake_val := sum(value), by= .(method, <<<redacted>>>, age_start, sex_id, year_id, draw)]
edu.data.long[variable %in% paste0("prop_",0:18), value := value /rake_val]
edu.data.long[, rake_val := NULL]

#insert means and weighted means from forecasted values
edu.data.long <- edu.data.long[variable != "weighted_mean"]
edu.data.long <- edu.data.long[!(variable == "mean" & ((year_id > year_before_forecasts & age_start <= 25)|
                                                         (year_id > year_before_forecasts + 5 & age_start == 30)|
                                                         (year_id>year_before_forecasts + 10 & age_start == 35)|
                                                         (year_id>year_before_forecasts + 15 & age_start == 40)|
                                                         (year_id>year_before_forecasts + 20& age_start == 45)|
                                                         (year_id>year_before_forecasts + 25& age_start == 50)|
                                                         (year_id>year_before_forecasts + 30& age_start == 55)|
                                                         (year_id>year_before_forecasts + 35& age_start == 60)|
                                                         (year_id>year_before_forecasts + 40& age_start == 65)|
                                                         (year_id>year_before_forecasts + 45& age_start == 70)|
                                                         (year_id>year_before_forecasts + 50& age_start == 75)|
                                                         (year_id>year_before_forecasts + 55& age_start == 80)|
                                                         (year_id>year_before_forecasts + 60& age_start == 85)|
                                                         (year_id>year_before_forecasts + 65& age_start == 90)|
                                                         (year_id>year_before_forecasts + 70& age_start == 95)))]


edu.data.long[, edu_yrs := as.numeric(str_sub(variable,6,-1))]
edu.wtdmean <- edu.data.long[variable %in% paste0("prop_",0:18),.(value = weighted.mean(edu_yrs, value)), by = .(method, <<<redacted>>>, age_start, sex_id, year_id,draw) ]
edu.wtdmean[, variable := "weighted_mean"]
edu.data.long<- rbind(edu.data.long, edu.wtdmean, fill = T)

edu.wtdmean[, variable := "mean"]
edu.wtdmean <- edu.wtdmean[(year_id > year_before_forecasts & age_start <= 25)|
                             (year_id > year_before_forecasts + 5 & age_start == 30)|
                             (year_id>year_before_forecasts + 10 & age_start == 35)|
                             (year_id>year_before_forecasts + 15 & age_start == 40)|
                             (year_id>year_before_forecasts + 20& age_start == 45)|
                             (year_id>year_before_forecasts + 25& age_start == 50)|
                             (year_id>year_before_forecasts + 30& age_start == 55)|
                             (year_id>year_before_forecasts + 35& age_start == 60)|
                             (year_id>year_before_forecasts + 40& age_start == 65)|
                             (year_id>year_before_forecasts + 45& age_start == 70)|
                             (year_id>year_before_forecasts + 50& age_start == 75)|
                             (year_id>year_before_forecasts + 55& age_start == 80)|
                             (year_id>year_before_forecasts + 60& age_start == 85)|
                             (year_id>year_before_forecasts + 65& age_start == 90)|
                             (year_id>year_before_forecasts + 70& age_start == 95)]
edu.data.long<- rbind(edu.data.long, edu.wtdmean, fill = T)


edu.data.long[,mean := NULL] %>%
  .[,edu_yrs2 := NULL] %>% 
  .[,logit_value := NULL]


#only keep final, single-year bins, means, and metrics of inequality, then reconstruct them from single year bins
edu.data.long <- edu.data.long[!variable%like% "cprop" & !variable %in% c("prop_1_5", "prop_6_11", "prop_12_18")]

#small numbers correction
edu.data.long[value <= .0005, value := .0005]

##make sure all variables are there
fin.data.wide <- dcast.data.table(data=edu.data.long,formula=<<<redacted>>>+age_start+sex_id+year_id+method+draw~variable,value.var = "value")
fin.data.wide[, prop_1_5 := sum(prop_1, prop_2, prop_3, prop_4, prop_5), by = .(<<<redacted>>>, age_start, sex_id, year_id, method, draw)]
fin.data.wide[, prop_6_11 := sum(prop_6, prop_7, prop_8, prop_9, prop_10, prop_11), by = .(<<<redacted>>>, age_start, sex_id, year_id, method, draw)]
fin.data.wide[, prop_12_18 := sum(prop_12, prop_13, prop_14, prop_15, prop_16, prop_17, prop_18), by = .(<<<redacted>>>, age_start, sex_id, year_id, method, draw)]

fin.data.wide[,cprop_0:=prop_0]
fin.data.wide[,cprop_12_18:=prop_12_18 / (1-cprop_0)]
fin.data.wide[,cprop_1_5:=prop_1_5 / (1-cprop_0+prop_12_18)]
#1-5
fin.data.wide[,cprop_1:=prop_1 / prop_1_5]
fin.data.wide[,cprop_2:=prop_2 / (prop_1_5-prop_1)]
fin.data.wide[,cprop_3:=prop_3 / (prop_1_5-(prop_1+prop_2))]
fin.data.wide[,cprop_4:=prop_4 / (prop_1_5-(prop_1+prop_2+prop_3))]
#6-11
fin.data.wide[,cprop_6:=prop_6 /prop_6_11]
fin.data.wide[,cprop_7:=prop_7 / (prop_6_11-prop_6)]
fin.data.wide[,cprop_8:=prop_8 / (prop_6_11-(prop_6+prop_7))]
fin.data.wide[,cprop_9:=prop_9 / (prop_6_11-(prop_6+prop_7+prop_8))]
fin.data.wide[,cprop_10:=prop_10 / (prop_6_11-(prop_6+prop_7+prop_8+prop_9))]
#12-18
fin.data.wide[,cprop_12:=prop_12 / prop_12_18]
fin.data.wide[,cprop_13:=prop_13 / (prop_12_18-prop_12)]
fin.data.wide[,cprop_14:=prop_14 / (prop_12_18-(prop_12+prop_13))]
fin.data.wide[,cprop_15:=prop_15 / (prop_12_18-(prop_12+prop_13+prop_14))]
fin.data.wide[,cprop_16:=prop_16 / (prop_12_18-(prop_12+prop_13+prop_14+prop_15))]
fin.data.wide[,cprop_17:=prop_17 / (prop_12_18-(prop_12+prop_13+prop_14+prop_15+prop_16))]

#reshape and merge on variable codes
fin.data.long <- melt.data.table(fin.data.wide, id.vars=c("method","<<<redacted>>>","age_start","sex_id","year_id","draw"))
fin.data.long[, code := variable]
edu.codes <- fread(paste0(root,"ref/edu_codes.csv"))
fin.data.long <- merge(fin.data.long,edu.codes,by=c("code"), all.x = T)


fin.data.long[is.na(variable), variable := code]
fin.data.long[is.na(code), code := variable]

#add aaid and gini variables
fin.data.long <- fin.data.long[code != "aaid"]
fin.data.long[, edu_yrs2 := as.numeric(str_sub(variable, 6, -1))]
means <- fin.data.long[code %in% paste0("prop_", 0:18),.(mean = weighted.mean(edu_yrs2, value)), by = .(method, <<<redacted>>>, age_start, sex_id, year_id,draw)]
fin.data.long <- merge(fin.data.long, means, by = c("method", "<<<redacted>>>", "age_start", "sex_id", 'year_id', "draw"), all.x = T)
aaid <- fin.data.long[variable %in% paste0("prop_", 0:18),.(aaid = ((gini(x = edu_yrs2, weights = value)) * 2 * mean)), by = .(draw,method, <<<redacted>>>, age_start, sex_id, year_id)]
aaid <- aaid[duplicated(aaid)==F]
aaid[,variable := "aaid"]
setnames(aaid, "aaid", "value")
fin.data.long <- rbind(aaid, fin.data.long, fill = T)
fin.data.long <- fin.data.long[,.(draw, method, <<<redacted>>>, age_start, sex_id, year_id, value, variable, edu_yrs, description, code)]

#reshape one final time, compute confidence intervals if working with draws, assert all quantities are there and save
edu.data.no.draws <- fin.data.long[,.(value = mean(value, na.rm = T), upper = quantile(value, .95, na.rm = T), lower = quantile(value, .05, na.rm = T)), by = .(method, <<<redacted>>>, age_start, sex_id, year_id, variable, edu_yrs, code, description)]

assert_ids(edu.data.no.draws, id_vars = list(age_start = seq(5,95,5),
                                             year_id = 1950:end_year,
                                             variable = c(paste0("prop_", c(0:18, "1_5", "6_11", "12_18")),
                                                          paste0("cprop_", c(0:4, 6:10, 12:17, "1_5", "12_18")),
                                                          "mean", "weighted_mean", "aaid")),
           assert_dups = T,
           warn_only = T)



#compress draws file
fin.data.long[,c("description", "method", "variable", "edu_yrs") := NULL]
fin.data.long[,code := as.factor(code)]
fin.data.long[,<<<redacted>>> := as.factor(<<<redacted>>>)]
fin.data.long[,draw := as.factor(draw)]
fin.data.long[,age_start := as.factor(age_start)]
fin.data.long[,sex_id := as.factor(sex_id)]
fin.data.long[,draw := as.factor(draw)]


wide <- dcast.data.table(fin.data.long, code + <<<redacted>>> + age_start + sex_id + year_id ~ draw, value.var = "value")
head(wide)
wide[code == "mean" & age_start == 75 & year_id == 1950]
#save data 
dir.create(paste0(hroot,"/data/output_data/",model_version,"/processed/",h.scen,"/"),recursive = T,showWarnings = F)
saveRDS(edu.data.no.draws,paste0(hroot,"/data/output_data/",model_version,"/processed/",h.scen,"/",c.loc_sex,".rds"))

#save processed data by code for fast raking
for(c.code in unique(edu.data.no.draws$code)){
  print(c.code)
  dir.create(paste0(hroot,"/data/output_data/",model_version,"/processed/",h.scen,"/", c.code),recursive = T,showWarnings = F)
  saveRDS(edu.data.no.draws[code == c.code],paste0(hroot,"/data/output_data/",model_version,"/processed/",h.scen,"/",c.code, "/",c.loc_sex,".rds"))
}
dir.create(paste0(hroot,"/data/output_data/",model_version,"/processed_draws/",h.scen,"/"),recursive = T,showWarnings = F)
saveRDS(wide,paste0(hroot,"/data/output_data/",model_version,"/processed_draws/",h.scen,"/",c.loc_sex,".rds"))
for(c.code in unique(wide$code)){
  print(c.code)
  dir.create(paste0(hroot,"/data/output_data/",model_version,"/processed_draws/",h.scen,"/", c.code),recursive = T,showWarnings = F)
  saveRDS(wide[code == c.code],paste0(hroot,"/data/output_data/",model_version,"/processed_draws/",h.scen,"/",c.code, "/",c.loc_sex,".rds"))
}
