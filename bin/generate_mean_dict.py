import sys

#  Import numpy for numerics
import numpy.ma as ma

import numpy as np

#  Import pylab for plotting
#import pylab as plt

#  Import Network class from pygeonetwork
#from pygeonetwork import Network

#  Import DynamicSocietyModel class from pydynamicsociety
#from pydynamicsociety import DynamicSocietyModel

# Import Networkx as a plotting tool for 
#import networkx as nx

#from pygeonetwork import mpi

import glob 
import os
import pickle
import scipy.stats


path =''#'/scratch/01/carls/social/emulation/ens_rc_p25_ens_100/'

#os.chdir(path)

run_reps=glob.glob('trans_smok*')
ensemble_size=len(run_reps)
ens_dic={}
full_dict={}
dict_disconnected={}

pkl_file=open(run_reps[0],'rb')            
read_in=pickle.load(pkl_file)
pkl_file.close()
fin_dic={}


for active_var in read_in[read_in.keys()[0]].iterkeys():
	fin_dic[active_var]={}

#for folder in run_reps:
#try:
# print folder
summary_dic={}
no_connected_list=[]
folder='ens'
dict_disconnected[folder]={}
for run in xrange(ensemble_size):
	pkl_file=open('trans_smok_%s.pkl'%(run),'rb')            
	read_in=pickle.load(pkl_file)
	print run
	for active_key in read_in.iterkeys():

		for active_var in read_in[active_key].iterkeys():
			#print active_var, active_key
			if active_var=='centrality_smokers_connected': 
				if active_key=='full':
					rel=read_in['full']['centrality']['smoker'][:10].mean()
					no_connected_list.append(read_in['full']['centrality_smokers_connected']/rel)	
				
			elif active_var=='disconnected_nodes':
				if active_key=='full':
					dict_disconnected[folder][run]=read_in[active_key][active_var]
				
			elif summary_dic.has_key(active_var) == False: 
				summary_dic[active_var]={}
			
			elif active_var=='smokers':
				#print 'smokers unpacked'
				if summary_dic[active_var].has_key(active_key) == False:
					summary_dic[active_var][active_key]=[]
				summary_dic[active_var][active_key].append(read_in[active_key][active_var])	
			
			elif np.asarray(read_in[active_key][active_var]).dtype=='object':
				for sm_key in read_in[active_key][active_var].iterkeys():
					if summary_dic[active_var].has_key(sm_key) == False:
						summary_dic[active_var][sm_key]={}
		
					if	summary_dic[active_var][sm_key].has_key(active_key)== False:
						summary_dic[active_var][sm_key][active_key]=[]
			
					summary_dic[active_var][sm_key][active_key].append(read_in[active_key][active_var][sm_key])
			else:	
				if summary_dic[active_var].has_key(active_key) == False:
					summary_dic[active_var][active_key]=[]
	
				summary_dic[active_var][active_key].append(read_in[active_key][active_var])	

			# if active_var=='centrality_smokers_connected' and active_key=='full': 
			# 	rel=read_in['full']['centrality']['smoker'][:10].mean()
			# 	no_connected_list.append(read_in['full']['centrality_smokers_connected'][-10:].mean()/rel)	
			

	pkl_file.close()			
	#centr_connected.append((folder,np.asarray(no_connected_list).mean(),np.asarray(no_connected_list).std()))
				
#############################
# generate mean
#############################
print 'generate mean'

ave_dic={}
for active_var in summary_dic.iterkeys():
	if ave_dic.has_key(active_var) == False: 
		ave_dic[active_var]={}
	
	#print active_var, active_key
	if active_var=='no_disconnected_nodes':
		for run_key in summary_dic[active_var].iterkeys():
			if ave_dic[active_var].has_key(run_key) == False: 
				ave_dic[active_var][run_key]={}
			dataset=np.asarray(summary_dic[active_var][run_key])
			dataset[dataset>100]=1
			#print active_var,dataset.shape
			ave_dic[active_var][run_key]['mean']=dataset.mean(axis=0)
			ave_dic[active_var][run_key]['lower']=dataset.min(axis=0)
			ave_dic[active_var][run_key]['upper']=dataset.max(axis=0)
				
	
	
	elif summary_dic[active_var].has_key('smoker'):
		#print 'has_key',active_var
		print active_var
		for active_key in summary_dic[active_var].iterkeys():
			if ave_dic[active_var].has_key(active_key) == False: 
				ave_dic[active_var][active_key]={}
		
			for run_key in summary_dic[active_var][active_key].iterkeys():
				if ave_dic[active_var][active_key].has_key(run_key) == False: 
					ave_dic[active_var][active_key][run_key]={}
					dataset=np.asarray(summary_dic[active_var][active_key][run_key])
					ave_dic[active_var][active_key][run_key]['mean']=scipy.stats.scoreatpercentile(dataset[:,:]/np.mean(dataset[:,1:25],axis=1).reshape(dataset.shape[0],1),50)
					ave_dic[active_var][active_key][run_key]['lower']=scipy.stats.scoreatpercentile(dataset[:,:]/np.mean(dataset[:,1:25],axis=1).reshape(dataset.shape[0],1),25)
					ave_dic[active_var][active_key][run_key]['upper']=scipy.stats.scoreatpercentile(dataset[:,:]/np.mean(dataset[:,1:25],axis=1).reshape(dataset.shape[0],1),75)
	
	elif np.asarray(summary_dic[active_var]['full']).dtype==object:
		print 'smokers yes'
		
		for run_key in summary_dic[active_var].iterkeys():
			if ave_dic[active_var].has_key(run_key) == False: 
				ave_dic[active_var][run_key]={}
			dataset=np.zeros((len(summary_dic[active_var][run_key]),len(summary_dic[active_var][run_key][0])))
			print "dataset.shape",dataset.shape
			for i in xrange(len(summary_dic[active_var][run_key])):
				for k in xrange(len(summary_dic[active_var][run_key][i])):
					#print summary_dic[active_var][run_key][i][k][0].shape[0]
					dataset[i,k]=summary_dic[active_var][run_key][i][k][0].shape[0]
	
			ave_dic[active_var][run_key]['mean']=dataset.mean(axis=0)
			ave_dic[active_var][run_key]['lower']=dataset.min(axis=0)
			ave_dic[active_var][run_key]['upper']=dataset.max(axis=0)
			
	else:
		for run_key in summary_dic[active_var].iterkeys():
			if ave_dic[active_var].has_key(run_key) == False: 
				ave_dic[active_var][run_key]={}
			dataset=np.asarray(summary_dic[active_var][run_key])
			#print active_var,dataset.shape
			ave_dic[active_var][run_key]['mean']=scipy.stats.scoreatpercentile(dataset,50)
			ave_dic[active_var][run_key]['lower']=scipy.stats.scoreatpercentile(dataset,25)
			ave_dic[active_var][run_key]['upper']=scipy.stats.scoreatpercentile(dataset,75)
	
#############################
# calculate final values
#############################

#print runsumdict
run_sum_dict=summary_dic
for active_var in run_sum_dict.iterkeys():
	if run_sum_dict[active_var].has_key('smoker'):
		#print 'has_key',active_var
		print active_var
		for active_key in run_sum_dict[active_var].iterkeys():
				dataset=np.asarray(run_sum_dict[active_var]['smoker']['full'])
				fin_dic[active_var][str(folder)[-2:]]=np.mean(dataset[:,-200:],axis=1)/np.mean(dataset[:,1:25],axis=1)
	elif np.asarray(run_sum_dict[active_var]['full']).dtype==object:
		print 'smokers yes'
		for run_key in run_sum_dict[active_var].iterkeys():
			dataset=np.zeros((len(run_sum_dict[active_var][run_key]),len(run_sum_dict[active_var][run_key][0])))
			print "dataset.shape",dataset.shape
			for i in xrange(len(run_sum_dict[active_var][run_key])):
				for k in xrange(len(run_sum_dict[active_var][run_key][i])):
					#print run_sum_dict[active_var][run_key][i][k][0].shape[0]
					dataset[i,k]=run_sum_dict[active_var][run_key][i][k][0].shape[0]
			fin_dic[active_var][str(folder)[-2:]]=np.mean(dataset[:,-200:],axis=1)

	else:
		for run_key in run_sum_dict[active_var].iterkeys():
			if fin_dic[active_var].has_key(run_key) == False: 
				fin_dic[active_var][run_key]={}
			dataset=np.asarray(run_sum_dict[active_var][run_key])
			#print active_var,dataset.shape
			fin_dic[active_var][str(folder)[-2:]]=np.mean(dataset[:,-200:],axis=1)



output = open('median_trans_smok.pkl', 'wb')
pickle.dump(ave_dic, output)
output.close()			
ens_dic=ave_dic
full_dict=summary_dic
	# except:
	# 	print 'generation failed for',folder

# output = open('ens_disconnected_nodes.pkl', 'wb')
# pickle.dump(dict_disconnected, output)
# output.close()
output = open('full_trans_smok.pkl', 'wb')
pickle.dump(full_dict, output)
output.close()

output = open('ens_trans_smok.pkl', 'wb')
pickle.dump(ens_dic, output)
output.close()			

output = open('final_reduction.pkl', 'wb')
pickle.dump(fin_dic, output)
output.close()

#np.savetxt('centr_connected_nodes_final.npy',centr_connected)
				
