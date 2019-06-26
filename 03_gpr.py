import os
os.chdir("/<<<<< filepath redacted >>>>>/")
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gpr.gpr as gpr
#reload(gpr)

def logit(p):
	return np.log(p) - np.log(1 - p)

def inv_logit(p):
	return np.exp(p) / (1 + np.exp(p))

#pass arguments
model_version = str(sys.argv[1])
c_draws = int(sys.argv[2])
h_schem = int(sys.argv[3])

#array task id tels us location-sex combination to model when running in parallel
task_id = int(os.environ.get("SGE_TASK_ID")) - 1
parameters = pd.read_csv('/<<<<< filepath redacted >>>>>/gpr_parameters.csv')
loc_sex = parameters.loc_sex[task_id]
c_iso = loc_sex[0:-2]
c_sex = loc_sex[-1]

#read in location set to merge on region varaibles
locs_set = pd.read_csv('<<<<< filepath redacted >>>>>/ref/locs.csv')
region_id = int(locs_set.region_id[locs_set['<<<redacted>>>_loc_id'] == c_iso])
reg_sex = str(region_id) + "_" + str(c_sex)

#load data
data = pd.read_csv('/<<<<< filepath redacted >>>>>/data/output_data/{model_version}/linear/{h_schem}/{reg_sex}.csv'.format(model_version=model_version,reg_sex=reg_sex,h_schem=h_schem))
data= data.loc[(data.<<<redacted>>>loc_id == '{c_iso}'.format(c_iso=c_iso))]
data['year'] = data['year_id']

df_list = []
gpr_data = pd.DataFrame()

#run GPR seperately for every age-quantity
for age in pd.unique(data['age_group_id']):
	for edyr in pd.unique(data['edu_yrs']):
		iso_sex_age_data = data.loc[((data.age_group_id == age) & (data.edu_yrs == edyr)),:]
		amp = iso_sex_age_data['mad'].values[0] * 1.4826 
		gpr_out = gpr.fit_gpr(iso_sex_age_data,year_variable='year',amp=amp,scale=40,obs_variable='logit_value', obs_var_variable='logit_variance', mean_variable='logit_value_pred',diff_degree=2,draws=c_draws)
		df_list.append(gpr_out)
		


print('test')

		
gpr_data = pd.concat(df_list)
gpr_data.to_csv('/<<<<< filepath redacted >>>>>/data/output_data/{model_version}/gpr/{h_schem}/{loc_sex}.csv'.format(model_version=model_version,loc_sex=loc_sex,h_schem=h_schem))




