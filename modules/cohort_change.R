#################################################################
#model age-cohort changes
### Joseph Friedman - josephfriedman@mednet.ucla.edu
### Hunter York - hyowrk@uw.edu
##################################################################
#subset cohort data
edu.all.subset <- edu.all[sample_size>100 & age_start>=25  ]
edu.all.subset[value<.05,value:=.05]
edu.all.subset[value>.95,value:=.95]
edu.all.subset[,cohort:=year-age_start]
##collapse cohorts to five years
edu.all.subset[, cohort := floor(cohort/5)*5]




cohort.info <- edu.all.subset[,.(min_age=min(age_start),max_age=max(age_start), count = length(value)),by=.(cohort,iso3,sex,type,edu_yrs)]

##only keep those with 2 or more observations
cohort.info <- cohort.info[min_age != max_age]
#subset to only DHS/ MICS for certain models
if (c.cohort.mod %in% c("3","4","5","6","7","8")) {cohort.info <- cohort.info[type %in% c("IPUMS","DHS")]}
edu.cohorts <- data.table(merge(edu.all.subset,cohort.info,by=c("cohort","iso3","sex","type","edu_yrs")))

##################################
##################################
#remove dups
edu.cohorts <- edu.cohorts[,.(value = mean(value)), by = .(cohort, iso3, sex, type, edu_yrs, age_start,   min_age, max_age)]

#create cohort variables
edu.cohorts[,coh.id := paste0(iso3, sex, type, edu_yrs, cohort)]

get_cohort_values <- function(c.coh.id){
  temp <- edu.cohorts[coh.id == c.coh.id]
  template <- expand.grid(cohort = unique(temp$cohort),
                          iso3 = (unique(temp$iso3)),
                          sex = unique(temp$sex),
                          type = unique(temp$type),
                          edu_yrs = unique(temp$edu_yrs),
                          orig_age = unique(temp$age_start),
                          new_age = unique(temp$age_start))
  template <- data.table(template)
  template <- template[new_age > orig_age]
  template[, orig_value := temp[age_start == orig_age, value], by = orig_age]
  template[, new_value := temp[age_start == new_age, value], by = new_age]
  template[, logit_diff := logit(new_value) - logit(orig_value)]
  return(template)
}

edu.cohorts.combs <- mclapply(unique(edu.cohorts$coh.id), get_cohort_values, mc.cores = cores)

edu.cohorts.diffs <- rbindlist(edu.cohorts.combs)
edu.cohorts.diffs[, age_start := new_age]
edu.cohorts.diffs <- merge(edu.all.subset, edu.cohorts.diffs, by = c("cohort", "iso3", "sex", "type", "edu_yrs", "age_start"))
edu.cohorts.diffs$new_value <- NULL
edu.cohorts.diffs[, age_diff := new_age - orig_age]
edu.cohorts.diffs$new_age <- NULL
edu.cohorts.diffs$orig_age <- NULL
edu.cohorts.diffs[, min_age := min(age_start), by = .(cohort, iso3, sex, type, edu_yrs)]
edu.cohorts.diffs[, max_age := max(age_start), by = .(cohort, iso3, sex, type, edu_yrs)]
locs <- fread(paste0(root,"<<<<< filepath redacted >>>>>/GBD_2019_locs_20181008.csv"))
locs2 <- locs[,c("region_name","super_region_name","<<<redacted>>>"),with=F]
locs2[,iso3:=<<<redacted>>>]
edu.cohorts.diffs <- data.table(merge(edu.cohorts.diffs,locs2,by=c("iso3")))
edu.cohorts.diffs[, sub_ihme_loc := substr(iso3,1,3)]
edu.cohorts <- edu.cohorts.diffs

#make prediction template
back_cast <- function(dt,time) {
  timecast.dt <- dt[age_start>=25+time]
  timecast.dt[,age_orig := age_start]
  timecast.dt <- timecast.dt[,age_start:= age_start - time]
  timecast.dt <- timecast.dt[,year:=year - time]
  timecast.dt[,age_diff:= -1*time]
  setnames(timecast.dt,"value","orig_value")
  #timecast.dt[,mean_se:=mean_se * ((time/5)+1)]
  return(timecast.dt)
  }

for_cast <- function(dt,time) {
  timecast.dt <- dt[age_start<=95-time & age_start>=25]
  timecast.dt[,age_orig := age_start]
  timecast.dt <- timecast.dt[,age_start:= age_start + time]
  timecast.dt <- timecast.dt[,year:=year + time]
  timecast.dt[,age_diff:=time]
  setnames(timecast.dt,"value","orig_value")
  #timecast.dt[,mean_se:=mean_se * ((time/5)+1)]
  return(timecast.dt)
  }

#Create Backcasts and Forecasts
back.casts.list <- mclapply(X=seq(5,55,5),FUN=back_cast,dt=edu.all,mc.cores=ifelse(Sys.info()[1]=="Windows", 1, cores))
for.casts.list <- mclapply(X=seq(5,70,5),FUN=for_cast,dt=edu.all,mc.cores=ifelse(Sys.info()[1]=="Windows", 1, cores))
template <- rbind(rbindlist(back.casts.list),rbindlist(for.casts.list))
template <- data.table(merge(template,locs2,by=c("iso3")))

#ensure edu.cohorts has age_orig variable
edu.cohorts[,age_orig := age_start-age_diff]

#Define Model and Predict Function



if(c.cohort.mod%in% c("10")) {
  ch.reg <-function(i) {
    
    print(i)
    c.edu.cohorts <- edu.cohorts[age_diff %in% 5:10 & year %in% 1990:2019]
    c.edu.cohorts[,age_mid_point := (age_orig + age_start)/2]
    c.edu.cohorts <- c.edu.cohorts[nchar(iso3)==3]
    c.edu.cohorts <- c.edu.cohorts[orig_value > .1 & value > .1]
    c.edu.cohorts <- c.edu.cohorts[orig_value < .9 & value < .9]
    c.edu.cohorts[,avg_diff := mean(logit_diff[age_mid_point %in% seq(60,70,.5)]), by = .(iso3, sex, type, edu_yrs, year, age_diff)]
    c.edu.cohorts[is.nan(avg_diff), avg_diff := 0]
    c.edu.cohorts[,logit_diff := logit_diff - avg_diff]
    c.edu.cohorts[,annualized_logit_diff := logit_diff/age_diff]
    mod <- lmer(annualized_logit_diff~ ns(age_mid_point,knots=c(75)) + (1|region_name/iso3),data=c.edu.cohorts[edu_yrs == i])
    temp <- template[edu_yrs == i]
    grid <- expand.grid(region_name = unique(c.edu.cohorts$region_name),iso3 = "XXX", age_mid_point = seq(25,95,1))
    grid <- data.table(grid)
    grid[, pred_diff := predict(mod,newdata=grid,allow.new.levels=T, re.form = ~(1|region_name))]
    grid[, confidint := (data.table(predictInterval(mod,
                                                    newdata=grid,
                                                    which = "fixed",
                                                    ignore.fixed.terms = c(1),
                                                    include.resid.var = F,
                                                    type = "linear.prediction")) %>% .[,(fit - lwr)])]
    
    #make sure nothing happens below 65
    grid[age_mid_point < 65, pred_diff := 0]
    
    expected <- data.table(expand.grid(region_name = unique(c.edu.cohorts$region_name), age_start = seq(25,95,5), age_end = seq(25,95,5)))
    expected <- expected[age_start != age_end]
    
    get_delta <- function(c.age_start, c.age_end, c.region_name){
      subset <- grid[age_mid_point %in% c.age_start:c.age_end & region_name == c.region_name]
      change <- sum(subset[, pred_diff])
      if(c.age_start > c.age_end){change <- change * -1}
      return(change)
    }
    get_ci <- function(c.age_start, c.age_end, c.region_name){
      subset <- grid[age_mid_point %in% c.age_start:c.age_end & region_name == c.region_name]
      ci <- sum(subset[, confidint])
      return(ci)
    }
    expected[, tot := get_delta(age_start, age_end, region_name), by =.(age_start, age_end, region_name) ]
    expected[, extrap_ci := get_ci(age_start, age_end, region_name), by =.(age_start, age_end, region_name) ]
    expected[, age_diff := age_end - age_start]
    expected[, age_orig := age_start]
    expected[,age_end := NULL]
    expected[, age_start := NULL]
    temp <- merge(temp, expected, by = c("region_name", "age_orig", "age_diff"), all.x = T)
    temp[, logit_diff_pred := tot] %>% .[,tot := NULL]
    return(temp) 
  } 
}

#Run Model and Predict in Parallel
ch.list <- lapply(unique(edu.cohorts$edu_yrs),ch.reg)
print(ch.list)
ch.pred<- rbindlist(ch.list,fill=T)

#Calculate predicted values
ch.pred[,value:=inv.logit(logit(orig_value) + logit_diff_pred)]

#Save Plots of Cohort Changes
ch.pred.small <- copy(ch.pred)

ch.pred.small[,change:=round(value-orig_value,2)]
ch.pred.small[,logit_diff_pred:=round(logit_diff_pred,2)]
ch.pred.small[,obs:=1]
ch.pred.small.tots <- ch.pred.small[,.(num=sum(obs)),by=.(super_region_name,edu_yrs)]
ch.pred.small <- ch.pred.small[,.(num=sum(obs)),by=.(change,logit_diff_pred,super_region_name,edu_yrs)]
ch.pred.small <- merge(ch.pred.small,ch.pred.small.tots,by=c("super_region_name","edu_yrs"))
ch.pred.small[,weight:=num.x/num.y]

ch.pred.small <- merge(ch.pred.small,edu.codes,by=c("edu_yrs"))

c.descs <- c("average years of schooling","proportion with 0","proportion of 1-18 with 12-18","proportion of 1-5 in 1"  )

c.desc <- "average years of schooling"

##create one dataset
edu.all[,point_type:=0]
ch.pred[, point_type := age_diff]

vars <-names(edu.all)
if(!"extrap_ci" %in% names(ch.pred)){
  ch.pred[, extrap_ci := 0]
}
ch.pred.graph <-  rbind(edu.all,ch.pred,fill=T)

ch.pred <- ch.pred[,c("iso3","year","age_start","sex","type","edu_yrs","sample_size","value","part1","part2","part3","part4","point_type","extrap_ci")]
edu.all.imp <- rbind(edu.all,ch.pred,fill=T)
edu.all.imp[is.na(extrap_ci), extrap_ci := 0]
