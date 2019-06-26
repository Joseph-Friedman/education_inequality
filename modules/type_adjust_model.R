### Joseph Friedman - josephfriedman@mednet.ucla.edu
### Hunter York - hyowrk@uw.edu
#prep data
locs2 <- fread(paste0(root,"ref/locs.csv"))
locs2[,iso3:=ihme_loc_id]
locs2<- locs2[,.SD,.SDcols=c("iso3","region_name","super_region_name","developed")]

#classify surveys
for(c.type in c("DHS","IPUMS","CENSUS","BRFSS","WHS","WFS","ISSP","EUROBAROMETER","EUROSTAT","MICS","RHS","MCSS","AFROBAROMETER","UKGHS", "LATINOBAROMETER", "ARAB_BAROMETER", "TRICT_LEVEL", "NSSO")){
  edu.all.imp[type %like% c.type, type2 := c.type]
}
edu.all.imp[, orig_year := year - point_type]
size <- edu.all.imp[point_type == 0,.(size=sum(sample_size,na.rm=T)),by=.(iso3,type,orig_year,edu_yrs)]
edu.all.imp <- merge(edu.all.imp,size,by=c("iso3","type","orig_year","edu_yrs"), all.x = T)
edu.all.imp[is.na(type2) & size < 10000, type2 := "other_small"]
edu.all.imp[is.na(type2) & size >= 10000, type2 := "other_large"]

# edu.all.imp[,type2:=factor(type2,levels=c("IPUMS",unique(edu.all.imp$type2)[unique(edu.all.imp$type2)!="IPUMS"]))]
edu.all.imp[,ts:=paste0(iso3,"_",age_start,"_",sex,"_",edu_yrs)]
edu.all.imp[,decade:=floor(year/10)*10]
edu.all.imp <- merge(edu.all.imp,locs2,by='iso3')

edu.all.imp[value == 0, value := .01]
edu.all.imp[value == 1, value := .99]
edu.all.imp[, logit_value := logit(value)]

###import country-specific gold standards 
country_gold_standards <- fread("<<<<< filepath redacted >>>>>/ref/gold_standards.csv", header = T, na.strings = "")
###
gt <- edu.all.imp
gt[, id := paste0(region_name, edu_yrs)]
adjust_data <- function(c.id){
  print(c.id)
  dat.2 <- copy(gt[id == c.id])
  dat.2[,sub_ihme := substr(iso3,1,3)]
  dat.2[,type_age := paste(type, age_start, sep = "_")]
  dat.2[,iso3_age := paste(iso3, age_start, sep = "_")]
  dat.2[,sub_type := substr(type2, 1,7)]
  dat.2[,type_iso3 := paste(sub_type, iso3, sep = "_")]
  dat.2[,standard_year := year-1990]
  dat.2[,standard_age := age_start-40]
  
  #######assign gold standard for region, either dhs or IPUMS
  type_counts <- dat.2[nchar(iso3)==3,.(count = length(unique(paste0(type2,sub_type,orig_year)))), by = type2]
  
  
  if(!"DHS" %in% unique(type_counts$type2)){
    type_counts <- rbind(type_counts,( data.table(table(type2 = "DHS", count = 0)) %>%  .[,N := NULL]))
  }
  if(!"IPUMS" %in% unique(type_counts$type2)){
    type_counts <- rbind(type_counts, data.table(table(type2 = "IPUMS", count = 0)) %>%  .[,N := NULL])
  }
  dat <- merge(dat.2,type_counts,by='type2')
  if(type_counts[type2 == "DHS", count] >=  
     type_counts[type2 == "IPUMS", count] & type_counts[type2 == "DHS", count] > 0){
    
    gold_standard <- "DHS"}else if(type_counts[type2 == "DHS", count] <  
        type_counts[type2 == "IPUMS", count] & type_counts[type2 == "IPUMS", count] > 0){
      gold_standard <- "IPUMS" 
    }else{
      gold_standard <- "NONE"
    }
  
  dat.2[,type2:=factor(type2,levels=c(gold_standard,unique(dat.2$type2)[unique(dat.2$type2)!=gold_standard]))]
  
  mod <- lmer(logit_value ~ iso3 +(1|sub_type/iso3) +standard_age+ standard_year+as.factor(sex):iso3, data = dat.2[age_start > 25 & age_start < 65 & nchar(iso3)==3])
  
  ranefs <- data.table(ranef(mod)$sub_type)
  ranefs[, sub_type:= rownames(ranef(mod)$sub_type)]
  ranefs[, logit_value := `(Intercept)`][, `(Intercept)`:= NULL]
  ipums_vals <- ranefs[sub_type == gold_standard]
  setnames(ipums_vals, "logit_value", "ipums_val")
  ipums_vals[, sub_type := NULL]
  ranefs[, ipums_val := ipums_vals]
  ranefs[, shift := ipums_val - logit_value]
  ranefs[is.na(shift), shift := 0]
  ranefs[, c("logit_value", "ipums_val"):=NULL]
  ranefs[, super_region_name := dat.2$super_region_name[1]]#[, sex := dat.2$sex[1]]
  dat.2 <- merge(dat.2, ranefs, by = c("sub_type", "super_region_name"), all.x = T)
  dat.2[,orig_value := logit_value]
  dat.2[, step_1_adj_logit_value := logit_value + shift]
  
  #add iso3-age shift
  ranefs <- data.table(ranef(mod)$`iso3:sub_type`)
  ranefs[,`iso3:sub_type`:= rownames(ranef(mod)$`iso3:sub_type`)]
  ranefs[, iso3 :=  tstrsplit(tstrsplit(`iso3:sub_type`, ":")[[1]], "_")[[1]]]
  ranefs[, sub_type := tstrsplit(`iso3:sub_type`, ":")[[2]]]
  ranefs[, logit_value := (`(Intercept)`)][, `(Intercept)`:= NULL][,`iso3:sub_type` :=NULL]
  
  ranefs <- merge(ranefs, country_gold_standards, by = "iso3", all.x = T)
  ranefs[sub_type == gold_standard & is.na(new_gold_standard) & iso3 != "SDN", new_gold_standard := gold_standard]
  ipums_vals <- ranefs[sub_type == new_gold_standard,.(ipums_val = mean(logit_value)), by = .(iso3)]
  other_vals <-  ranefs[,.(logit_value = mean(logit_value)), by = .(iso3,sub_type)]
  all_vals <- merge(ipums_vals, other_vals, by = "iso3")
  all_vals[,diff := ipums_val-logit_value][,logit_value := NULL]
  #avgs_by_type <- all_vals[,.(avg_change = mean(diff)), by = .(sub_type)]
  ##find countries with ipums and merge onto them their actual values, for all other places use average
  ipums_locs <- unique(ipums_vals$iso3)
  #dat.2 <- merge(dat.2, avgs_by_type, by = c("sub_type"), all.x = T)
  dat.2 <- merge(dat.2, all_vals, by.x = c("sub_type", "sub_ihme"), by.y = c("sub_type", "iso3"), all.x = T)
  dat.2[substr(iso3,1,3) %in% ipums_locs, final_shift := diff]
  #dat.2[!iso3 %in% ipums_locs, final_shift := 0]
  dat.2[is.na(final_shift), final_shift := 0]
  
  dat.2[, logit_value := step_1_adj_logit_value + final_shift]

  return(dat.2)
}

##save gold standard
get_gold_standard <- function(c.id){
  print(c.id)
  dat.2 <- copy(gt[id == c.id])
  type_counts <- dat.2[nchar(iso3)==3,.(count = length(unique(paste0(type2,iso3,orig_year)))), by = type2]
  
  
  if(!"DHS" %in% unique(type_counts$type2)){
    type_counts <- rbind(type_counts,( data.table(table(type2 = "DHS", count = 0)) %>%  .[,N := NULL]))
  }
  if(!"IPUMS" %in% unique(type_counts$type2)){
    type_counts <- rbind(type_counts, data.table(table(type2 = "IPUMS", count = 0)) %>%  .[,N := NULL])
  }
  dat <- merge(dat.2,type_counts,by='type2')
  if(type_counts[type2 == "DHS", count] >=  
     type_counts[type2 == "IPUMS", count] & type_counts[type2 == "DHS", count] > 0){
    
    gold_standard <- "DHS"}else if(type_counts[type2 == "DHS", count] <  
                                   type_counts[type2 == "IPUMS", count] & type_counts[type2 == "IPUMS", count] > 0){
      gold_standard <- "IPUMS" 
    }else{
      gold_standard <- "NONE"
    }
  c.reg <- unique(dat.2$region_name)
  output<- data.table(gold_standard = gold_standard, region_name = c.reg)
  setnames(output, "gold_standard", "region_gold_standard")
  
  c.country_gs <- country_gold_standards[iso3 %in% unique(substr(dat.2$iso3,1,3))]
  c.country_gs[, region_gold_standard := output$region_gold_standard]
  c.country_gs[is.na(new_gold_standard), new_gold_standard := "NONE"]
  
  return(c.country_gs)
}
gold_standard_list <- lapply(unique(gt$id), get_gold_standard)
gold_standards<- rbindlist(gold_standard_list) %>% unique

dir.create(paste0(hroot,"/data/output_data/",model_version,"/gold_standards/"), recursive = T)
write.csv(gold_standards, paste0(hroot,"/data/output_data/",model_version,"/gold_standards/",h.scen,".csv"), row.names = F)

##run source adj
#adj_list <- mclapply(unique(gt$id), adjust_data, mc.cores = cores)
adj_list <- lapply(unique(gt$id), adjust_data)
adj_gt <- rbindlist(adj_list)
edu.all.imp <- copy(adj_gt)


edu.all.imp[,value := inv.logit(logit_value)]
