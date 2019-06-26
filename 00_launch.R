################################################
### Launch Education Code Cascade
### Joseph Friedman - josephfriedman@mednet.ucla.edu
### Hunter York - hyowrk@uw.edu
################################################
list.of.packages <- c("data.table","stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
lapply(list.of.packages, require, character.only = TRUE)
#Arguments
jpath <- <<<<< filepath redacted >>>>>
root <- paste0(<<<<< filepath redacted >>>>>)
code_root <- paste0(<<<<< filepath redacted >>>>>)
################################################
### Set Run Settings
################################################
c.proj <- "proj_covariates"
print(commandArgs())
model_version <- commandArgs()[3]
c.cohort.mod <- as.numeric(str_sub(commandArgs()[4],1,-1))
c.source.adj.mod <- as.numeric(str_sub(commandArgs()[5],1,-1))
forecasting <- str_sub(commandArgs()[6],1,-1)
model_version <- "model_name"
c.cohort.mod <- "10"
c.source.adj.mod <- "4"
c.end_year <- 2040
forecasting <- "none"
queue <- "all.q"
cluster <- "new"
#model version settings
c.dists <- 1
c.draws <- 0
c.pv <- 0
c.pv2 <- 0

#which parts of code cascade to run
r.data  <- 1
r.model <- 1
r.gpr   <- 1
r.process  <- 1


#holdouts iterator 
if (c.pv == 1) {
  h.scens <- c(seq(0,9,1), 13,14)
} else {
  h.scens <- 0
}
if (c.pv2 == 1) {
  h.scens <- c(h.scens, seq(10,12,1))
} 

################################################
###load region-sex groups
################################################
locs <- fread(paste0(root,".csv")) #location set information
if(c.subnats == 0){locs <- locs[level==3]}
regs <- unique(locs$region_id)
reg.sexs <- c(paste0(regs,"_1"),paste0(regs,"_2"))

########################################################
# 01 - Data Prep
#######################################################
if (r.data==1) {
  for (h.scen in h.scens) {
    system(paste0(<<<<< redacted >>>>>))
  }
}

########################################################
# 02 - Model
#######################################################
if (r.model==1) {
  for (h.scen in h.scens) {
    for (c.reg_sex in reg.sexs) {
      system(paste0(<<<<< redacted >>>>>)
    }
  }
}

########################################################
# 03 - GPR
#######################################################


####create array jobs parameters map
loc_sexs <- c(paste0(unique(locs[level >= 3, ihme_loc_id]), "_1"), 
              paste0(unique(locs[level >= 3, ihme_loc_id]), "_2"))
n_jobs <- length(loc_sexs)
print(n_jobs)
Sys.sleep(2)
## Save the parameters as a csv so then you can index the rows to find the appropriate parameters
param_map <- expand.grid(loc_sex = loc_sexs)
write.csv(param_map, paste0(code_root,"/gpr_parameters.csv"), row.names=F)

####submit jobs
if (r.gpr==1) {
  for (h.scen in h.scens) {
    system(paste0(<<<<< redacted >>>>>))
  }
}
########################################################
#######################################################

if (r.process==1) {
  for (h.scen in h.scens) {
      system(paste0(<<<<< redacted >>>>>))
  }
}


########################################################
# 99 - PV
#######################################################
if (r.pv==1) {
  for (h.scen in h.scens) {
    system(paste0(<<<<< redacted >>>>>)
  }
}

