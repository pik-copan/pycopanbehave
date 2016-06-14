import glob 
import os
import pickle
from rpy2.robjects.vectors import DataFrame, StrVector
import rpy2.robjects as ro
import numpy as np
import subprocess

output = open('ens_trans_smok.pkl', 'rb')
ens_dic=pickle.load(output)
output.close()

no_connected_list=[]
for run_key in ens_dic.iterkeys():
	print run_key
	ave_dic=ens_dic[run_key]
	subprocess.call('mkdir %s'%run_key, shell=True)
	var_key='no_disconnected_nodes'
	for var_key in ave_dic.iterkeys():
		print 'var_key',var_key
		if ave_dic[var_key].has_key('smoker'):
			steps=ave_dic[var_key]['smoker']['dyn']['mean'].shape[0]
			mean=np.concatenate((ave_dic[var_key]['smoker']['dyn']['mean'],ave_dic[var_key]['smoker']['infl']['mean'],ave_dic[var_key]['smoker']['full']['mean']))
			lower=np.concatenate((ave_dic[var_key]['smoker']['dyn']['lower'],ave_dic[var_key]['smoker']['infl']['lower'],ave_dic[var_key]['smoker']['full']['lower']))
			upper=np.concatenate((ave_dic[var_key]['smoker']['dyn']['upper'],ave_dic[var_key]['smoker']['infl']['upper'],ave_dic[var_key]['smoker']['full']['upper']))
			timestep=np.concatenate((np.arange(1,steps+1),np.arange(1,steps+1),np.arange(1,steps+1)))
			condition=np.concatenate((np.repeat('smoker_dyn',steps),np.repeat('smoker_infl',steps),np.repeat('smoker_full',steps)))
			
			df_plot=DataFrame({'timestep': ro.IntVector(timestep),'condition' : StrVector(condition), 'mean':ro.FloatVector(mean) , 'lower':ro.FloatVector(lower),'upper':ro.FloatVector(upper)})
			df_plot.to_csvfile('%s/%s.csv'%(run_key,var_key),row_names=False,sep=' ')

			########################################################################
			# GGPLOT2 ENVIRONMENT IS USED AS IN THE BASIS R SCRIPT AND PORTED TO RPY2
			########################################################################

			subprocess.call('Rscript ggplot2_var_smok_dep.R %s %s %s'%(var_key,run_key,"Relative"), shell=True)
		elif var_key=='centrality_smokers_connected': 
			print var_key
		else:
			steps=ave_dic[var_key]['dyn']['mean'].shape[0]
			mean=np.concatenate((ave_dic[var_key]['dyn']['mean'],ave_dic[var_key]['infl']['mean'],ave_dic[var_key]['full']['mean']))
			lower=np.concatenate((ave_dic[var_key]['dyn']['lower'],ave_dic[var_key]['infl']['lower'],ave_dic[var_key]['full']['lower']))
			upper=np.concatenate((ave_dic[var_key]['dyn']['upper'],ave_dic[var_key]['infl']['upper'],ave_dic[var_key]['full']['upper']))
			timestep=np.concatenate((np.arange(1,steps+1),np.arange(1,steps+1),np.arange(1,steps+1)))
			condition=np.concatenate((np.repeat('smoker_dyn',steps),np.repeat('smoker_infl',steps),np.repeat('smoker_full',steps)))


			df_plot=DataFrame({'timestep': ro.IntVector(timestep),'condition' : StrVector(condition), 'mean':ro.FloatVector(mean) , 'lower':ro.FloatVector(lower),'upper':ro.FloatVector(upper)})
			df_plot.to_csvfile('%s/%s.csv'%(run_key,var_key),row_names=False,sep=' ')

			########################################################################
			# GGPLOT2 ENVIRONMENT IS USED AS IN THE BASIS R SCRIPT AND PORTED TO RPY2
			########################################################################

			subprocess.call('Rscript ggplot2_var_smok_dep.R %s %s %s'%(var_key,run_key,"Absolute"), shell=True)
			