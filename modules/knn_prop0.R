##KNN Module; takes in modeled mean and prop 0 to create full distribution of education using K nearest neighbors algorithm
### Joseph Friedman - josephfriedman@mednet.ucla.edu
### Hunter York - hyowrk@uw.edu

library(reldist, lib.loc = "<<<<< filepath redacted >>>>>")

#set hyperparameters
training_size <- 80
space_weight <- .8
cohort_weight <- .9
age_weight <- .15
c.span <- .5
psi <- 2.5

library(RcppZiggurat, lib.loc = "<<<<< filepath redacted >>>>>")
if (sessionInfo()$R.version$version.string %like% "3.5.1") {library(Rfast, lib.loc = "<<<<< filepath redacted >>>>>")}
library(reldist, lib.loc = "<<<<< filepath redacted >>>>>")


if("Rfast" %in% installed.packages()){yes_rfast <- TRUE}else{yes_rfast <- FALSE}
print(yes_rfast)

#load training_set
all.train <-readRDS(paste0(hroot,"/data/output_data/",model_version,"/input_data/", h.scen,"/raw_input_data.rds"))
if("orig_value" %in% names(all.train)){all.train[, value := orig_value]}
all.train <- all.train[code %in% paste0("prop_", 0:18)]

all.train[,year_id := year]
all.train[,sex_id := sex]
all.train[,<<<redacted>>> := iso3]
all.train[,cohort := year_id - age_start]
all.train[,edu_yrs := as.numeric(str_sub(code, 6, -1))]
all.train <- all.train[sample_size > 250]
all.train<- all.train[,.(value = mean(value)), by = .(<<<redacted>>>, iso3, type, year, cohort, sex_id, age_start, sex, edu_yrs)]
all.train[, trainer_id := paste(<<<redacted>>>, type, year,cohort, sex_id, sep = "_")]
all.train[,mean := weighted.mean(edu_yrs, value, na.rm =T), by = .(iso3, year, age_start, sex,type, trainer_id)]
####

all.train <- all.train[age_start %in% 5:100]
locs <- fread("<<<<< filepath redacted >>>>>/locs.csv")
all.train[, sub_ihme_loc := substr(<<<redacted>>>,1,3)]

all.train[,code := paste0("prop_", edu_yrs)]

#get aaid for each survey-year-loc-sex
prop0 <- all.train
prop0 <- prop0[,.(prop0 = value[code == "prop_0"]),by=.(age_start, year, sex_id, type, <<<redacted>>>, trainer_id)] #weighted gini; x is bin number, w is proportion in bin
all.train<-merge(all.train,prop0, by = c("age_start", "year", "sex_id", "type", "<<<redacted>>>", "trainer_id"), all.x = T)


#clean up 0 and 1 values for logit space
all.train[prop0 == 0, prop0 := 0.001]
all.train <- all.train[mean != 18 & mean != 0]
all.train <- data.table(merge(all.train, locs, by = c("<<<redacted>>>")))
all.train[, cohort := year- age_start]

#####
all.train[, sub_ihme_loc := substr(<<<redacted>>>,1,3)]
all.train[, location_id := NULL]
all.train <- merge(locs[, c("<<<redacted>>>", "location_id"), with = F], all.train, by.y = "sub_ihme_loc", by.x = "<<<redacted>>>")
setnames(all.train, "year", "year_id")

all.train.table <- all.train[value > 0.0001 ,.(count = length(unique(edu_yrs)),
                                 max = max(edu_yrs)), by = .(trainer_id, age_start)]
all.train <- all.train[!(trainer_id %in% all.train.table[count < 6 & age_start >= 15 & max > 6]$trainer_id)]
all.train.table <- NULL


all.train<- all.train[, c("prop0", "type", "value","age_start", "year_id", "mean", "trainer_id", "edu_yrs", "cohort", "location_id", "sex_id", "region_id", "super_region_id"), with = F]


all.train[, trainer_id := as.factor(trainer_id)]
all.train[, trainer_id := as.numeric(trainer_id)]
setkey(all.train, "trainer_id")
head(all.train)


all.train[age_start == 5, mean := mean / 3]
all.train[age_start == 10, mean := mean / 8]
all.train[age_start == 15, mean := mean / 13]
all.train[age_start > 15, mean := mean / 18]


all.train[, logit_mean := logit(mean)]
all.train[, logit_prop0 := logit(prop0)]

all.train <- all.train[!is.na(logit_mean) & !is.na(logit_prop0)]




c.all.train <- all.train[,.(proportion = sum(value)),by =  .( edu_yrs, trainer_id, sex_id, cohort, location_id, age_start, prop0)]
c.all.train[, mean := weighted.mean(edu_yrs, proportion), by = .(trainer_id, sex_id)]

means <- c.all.train[,.(mean = mean(mean)), by = .(trainer_id, sex_id)]
template <- data.table(expand.grid(edu_yrs = 0:18, trainer_id = unique(c.all.train$trainer_id)))
template_means <- merge(means, template, by = c("trainer_id"), all.y = T)
c.all.train <- merge(c.all.train[, .(edu_yrs, trainer_id, proportion)], template_means, all.y = T, by = c("trainer_id", "edu_yrs"))
c.all.train[is.na(proportion), proportion := 0]

setkey(c.all.train, "trainer_id")
all.train.males<- c.all.train[sex_id == 1] %>% .[,sex_id := NULL]
all.train.females<- c.all.train[sex_id == 2] %>% .[,sex_id := NULL]

all.rank <- all.train[,.(logit_mean = logit(mean(mean)), logit_prop0 = logit(mean(prop0))), by = .(age_start, year_id, trainer_id, cohort, location_id, sex_id, region_id, super_region_id)]
setkey(all.rank, "trainer_id")
rank.males<- all.rank[sex_id == 1]
rank.females<- all.rank[sex_id == 2]



load_data <- function(data){
  #loop over country-year-draw
  data[,iso_year_age_sex:=paste0(<<<redacted>>>,"_",year_id, "_",age_start, "_", sex_id, "_", method, "_", draw)]
  return(data)
}
get_dists <- function(data) {  
  #identify training set
  print(data)
  cc.iso <- unique(data$location_id)
  target.r <- unique(locs[location_id==cc.iso,region_id])
  target.sr <- unique(locs[location_id==cc.iso,super_region_id])
  if(data[, sex_id] == 1){train <- all.train.males} else {train <- all.train.females}
  if(data[, sex_id] == 1){c.rank <- rank.males} else {c.rank <- rank.females}
  c.rank[,spatial_distance:=1]
  c.rank[location_id %in% unique(locs[super_region_id==target.sr,location_id]),spatial_distance:=.66]
  c.rank[location_id %in% unique(locs[region_id==target.r,location_id]),spatial_distance:=.33]
  c.rank[location_id== cc.iso,spatial_distance:=0]
  c.year  <- data$year_id
  c.rank[, age_distance := abs(age_start - data$age_start)/75]
  c.rank[,cohort_distance := abs(cohort - data$cohort)/137]
  if (data$age_start > 20 & data[,mean] > .16){c.rank <- c.rank[age_start > 20]}else if (data$age_start ==5){c.rank <- c.rank[age_start == 5]} else if (data$age_start ==10){c.rank <- c.rank[age_start == 10]}# else if (data$age_start ==15){c.rank <- c.rank[age_start == 15]} else if (data$age_start ==20 & data$mean > .27){c.rank <- c.rank[age_start == 20]}else if (data$age_start ==20 & data$mean <= .27){c.rank <- c.rank[age_start <= 20]}
  
  
  
  c.rank[, distance := (spatial_distance * space_weight)^psi + (cohort_distance * cohort_weight)^psi + (age_distance * age_weight)^psi]
  
  empir_var <- var(c.rank[, list(logit_mean, logit_prop0)])
  
  if(yes_rfast){
    mahalo_distance <- mahala(x = as.matrix(c.rank[, list(logit_mean, logit_prop0)]),
                              mu = as.matrix(data[, list(logit(mean), logit(prop_0))]),
                              sigma = empir_var)
  } else{
    mahalo_distance <- mahalanobis(x = as.matrix(c.rank[, list(logit_mean, logit_prop0)]),
                                   center = as.matrix(data[, list(logit(mean), logit(prop_0))]),
                                   cov = empir_var)
  }
  
  c.rank[, mahalo_distance := mahalo_distance]
  
  c.rank[,final_rank:=frank(mahalo_distance,ties.method = "random")]
  c.rank[, mahalo_distance := NULL]
  
  #only subset to places with mean ed near the data, define cutoff
  train <- train[trainer_id %in% c(unique(c.rank[final_rank %in% 1:training_size]$trainer_id),
                                   unique(c.rank[location_id == cc.iso & inv.logit(logit_mean) < (data$mean + .02) & inv.logit(logit_mean) > (data$mean - .02) &
                                                   inv.logit(logit_prop0) < (data$prop_0 + .02) & inv.logit(logit_prop0) > (data$prop_0 - .02)]$trainer_id))]
  
  ##get ready to collapse
  #create vector of surveys used
  #collapse
  inv_distance <- c.rank[trainer_id %in% unique(train$trainer_id), .(trainer_id, distance)]
  train <- merge(train, inv_distance)
  train[distance < .01, distance := .01]
  train[, inv_distance := (1/distance)]
  
  new_props <- train[,.(proportion = weighted.mean(proportion,inv_distance)), by = .(edu_yrs)]
  new_props[,year_id :=  data[,year_id]]
  new_props[,sex_id := data[,sex_id]]
  new_props[,age_start := data[,age_start]]
  new_props[,draw := data[,draw]]
  return(new_props)
}

smooth_and_rake <- function(proportions){
  #fit loess smoother for each line
  c.proportions <- copy(proportions)
  c.proportions[, edu_yrs_sex_age_start := paste(edu_yrs, sex_id, age_start,method,draw, sep ="_")]
  c.proportions[new_proportion == 0, new_proportion:= .0001]
  c.proportions[new_proportion == 1, new_proportion:= .9999]
  
  smoothr <- function(demographic){
    c.props <- copy(c.proportions[edu_yrs_sex_age_start == demographic])
    lo <-  loess(logit(new_proportion) ~ year_id, degree = 2,span = c.span, data = c.props)
    c.props[, smooth := inv.logit(predict(lo))]
    return(c.props)
  }
  
  c.proportions <- rbindlist(mclapply(unique(c.proportions$edu_yrs_sex_age_start), smoothr, mc.cores = cores))
  
  #force the bins back to 0 if necessary 
  c.proportions[age_start ==5 & edu_yrs > 3, smooth := 0]
  c.proportions[age_start ==10 & edu_yrs > 8, smooth := 0]
  c.proportions[age_start ==15 & edu_yrs > 13, smooth := 0]
  
  #rake after the loess smoother
  c.proportions[, second_rake_factor := sum(smooth), by = .(year_id, <<<redacted>>>, sex_id, age_start, method,draw)]
  c.proportions[, final_estimate:= smooth/second_rake_factor]
  
  #compare to actual estimates
  c.proportions[, weighted_mean := weighted.mean(edu_yrs, final_estimate), by = .(<<<redacted>>>, age_start, year_id, sex_id, method,draw)]
  
  return(c.proportions)
}


#
locs <- fread(paste0(root,"ref/locs.csv"))
edu.dat_prepped <- load_data(edu.dat)

edu.dat_prepped[, sub_ihme_loc := substr(<<<redacted>>>, 1,3)]
edu.dat_prepped[, cohort := year_id - age_start]
edu.dat_prepped <- edu.dat_prepped[method == "gpr_mean"]
edu.dat_prepped <- merge(locs, edu.dat_prepped, by.x = "<<<redacted>>>", by.y = "sub_ihme_loc")
edu.dat_prepped <- edu.dat_prepped[, c("location_id", "age_start", "sex_id", "year_id", "method", "prop_0", "mean", "iso_year_age_sex", "cohort", "<<<redacted>>>", "draw")]
edu.dat_prepped[, draw := as.numeric(str_sub(draw,5,-1))]
edu.dat_prepped[, <<<redacted>>> := NULL]
edu.dat_prepped[, method := NULL]

#subset data to list of means to match on
subsettr <- function(c.id){print(c.id)
  subsetted <- edu.dat_prepped[iso_year_age_sex == c.id]
  subsetted$iso_year_age_sex <- NULL
  return(subsetted)}

prepped_data_list <- mclapply(unique(edu.dat_prepped$iso_year_age_sex), subsettr, mc.cores = cores)

#run KNN on each item of list
knn_dists_rough <- rbindlist(mclapply(prepped_data_list,get_dists, mc.cores = cores))

#set restrictions for lower ages
knn_dists_rough[age_start ==5 & edu_yrs > 7, proportion := 0]
knn_dists_rough[age_start ==10 & edu_yrs > 12, proportion := 0]
knn_dists_rough[age_start ==15 & edu_yrs > 17, proportion := 0]

#
knn_dists_rough[, draw := paste0("draw", draw)]
knn_dists_rough[, <<<redacted>>> := unique(edu.dat$<<<redacted>>>)]
knn_dists_rough[, method := unique(edu.dat$method)]

#rake a first time
knn_dists_rough[, rake_factor := sum(proportion), by = .(year_id, sex_id, <<<redacted>>>, age_start, method,draw)]
knn_dists_rough[, new_proportion := proportion/rake_factor]

#loess smooth the distributions and rake again
knn_dists_smooth <- smooth_and_rake(knn_dists_rough)

#
setnames(knn_dists_smooth, "weighted_mean", "mean")
knn_dists <- knn_dists_smooth[, c("edu_yrs", "year_id", "sex_id", "<<<redacted>>>", "age_start", "final_estimate", "mean", "method","draw")]
knn_dists <- data.table(dcast(knn_dists,year_id + sex_id + <<<redacted>>> + age_start + method +draw+ mean ~ edu_yrs, value.var = c("final_estimate")))

knn_dists[,prop_1_5:=`1` + `2` +`3` +`4` +`5`]
knn_dists[,prop_6_11:=`6` + `7` +`8` +`9` +`10` +`11`]
knn_dists[,prop_12_18:=`12` + `13` +`14` +`15` +`16` +`17`+`18`]
#entities used to model larger bins
knn_dists[,cprop_0:=`0`]
knn_dists[,cprop_12_18:=prop_12_18 / (1-cprop_0)]
knn_dists[,cprop_1_5:=prop_1_5 / (1-cprop_0+prop_12_18)]
#calculate conditional probabilities to model single year bins
#1-5
knn_dists[,cprop_1:=`1` / prop_1_5]
knn_dists[,cprop_2:=`2` / (prop_1_5-`1`)]
knn_dists[,cprop_3:=`3` / (prop_1_5-(`1`+`2`))]
knn_dists[,cprop_4:=`4` / (prop_1_5-(`1`+`2`+`3`))]
#6-11
knn_dists[,cprop_6:=`6` /prop_6_11]
knn_dists[,cprop_7:=`7` / (prop_6_11-`6`)]
knn_dists[,cprop_8:=`8` / (prop_6_11-(`6`+`7`))]
knn_dists[,cprop_9:=`9` / (prop_6_11-(`6`+`7`+`8`))]
knn_dists[,cprop_10:=`10` / (prop_6_11-(`6`+`7`+`8`+`9`))]
#12-18
knn_dists[,cprop_12:=`12` / prop_12_18]
knn_dists[,cprop_13:=`13` / (prop_12_18-`12`)]
knn_dists[,cprop_14:=`14` / (prop_12_18-(`12`+`13`))]
knn_dists[,cprop_15:=`15` / (prop_12_18-(`12`+`13`+`14`))]
knn_dists[,cprop_16:=`16` / (prop_12_18-(`12`+`13`+`14`+`15`))]
knn_dists[,cprop_17:=`17` / (prop_12_18-(`12`+`13`+`14`+`15`+`16`))]

for(i in as.character(0:18)) {setnames(knn_dists, i, paste0("prop_", i))}

knn_dists[, dist_type := "knn"]
setnames(knn_dists, "mean", "weighted_mean")
knn_dists[, draw := as.factor(draw)]
edu.dat$cprop_0 <- NULL
edu.dat$prop_0 <- NULL

edu.dat <- merge(knn_dists, edu.dat, by = c("year_id", "sex_id", "<<<redacted>>>", "age_start", "method", "draw"), all.y = T)
