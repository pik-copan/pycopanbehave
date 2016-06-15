#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
# Smoking Transition -- A plug for the pynamic society model to emulate smoking behaviour
#
# Copyright (C) 2013 Potsdam Institute for Climate Impact Research
# Authors: Jonathan F. Donges <donges@pik-potsdam.de>,
#          Carl-Friedrich Schleussner <schleussner@pik-potsdam.de>
# https://github.com/pik-copan/pycopanbehave

MAIN FUNCTIONS FOR SMOKING BEHAVIOUR TRANSITION SIMULATIONS

"""

#
#  Imports
#

#  Import numpy for numerics
import numpy as np

#Import os for system processes
import os

# Import igraph for network toolbox
import igraph

#  Import progressbar
import progressbar

import pickle

import scipy.stats
import scipy.integrate as integ

# Add repository path of the model core
import sys

sys.path.append('../')

# Import pyunicorn dependencies
from pyunicorn import mpi,Network

#  Import DynamicSocietyModel class from pydynamicsociety
from pysoc import DynamicSocietyModel

class Bunch(dict):
    """Bunch object"""
    def __init__(self,**kw):
        dict.__init__(self,kw)
        self.__dict__.update(kw)

####################################################################################################################################
#
#  Functions
#
####################################################################################################################################

def distribution_function(yb,x):
	"""" modulate the parabolic distribution y= a(b-x)**2 +c fct. with the following boundary conditions
			1. x=0: y== 3 : follows from the initial distribution function
			2. Int[0,1]==1: normalization criteria
			3. x=1: y==yb, where yb is the tunable parameter, which reflects the probability for an agent to have smoking disposition 1 
			yb0==3
	"""
	# print 'yb',yb
	#yb=3
	b=(1.+yb/3)/(1+yb)
	a=2./(b-1./3)
	c=3.-a*b**2

	return   a*(b-x)**2 +c

def rejection_sampling(N, distribution_function_ini,yb):
	#  Creates random sampling of an arbitraty distribution using the rejection sampling method for the continous characteristic (eg. smoking affinity) and a second binary characteristic based on this
	result = np.zeros((2,N))
	i = 0
	while i < N:
		random_numbers = np.random.rand(2)
		if random_numbers[0] < distribution_function_ini(yb,random_numbers[1]):
			result[0,i] = random_numbers[1]
			# if the smoking preference is greater than a certain value, the binary characteristic smoker/non-smoker is assigned
			result[1,i]=(random_numbers[1] > .5).astype('int8')
			i+=1
	return result


def interaction_likelihood_function(distance_metric_matrix,L):
	"""
	Returns the interaction_likelyhood matrix given a distance metric on the
	acquaintance network.
	The interaction likelyhood is treated exponential with a exp(-d/9)-relationship according 
	to the three degrees of influence 

	Since the number of nodes with a certain path length is gaussian with a peak around the mean apl, 
	this has to be accounted for. Thereby, we norm the interaction likelyhood according to the 
	absolute number of occurences of this it's entry

	"""
	# three degrees of influence plus an interaction offset to avoid disconnection
	exp_dec=(L.p_ai-L.interaction_offset)*np.exp(-(distance_metric_matrix-1)/2.)+L.interaction_offset
	distmax=distance_metric_matrix[np.isfinite(distance_metric_matrix)].max()
	histo_range=np.arange(1,distmax)
	distribution=np.histogram(distance_metric_matrix.flatten(),histo_range)
	for i in distribution[1][:-1]:
		exp_dec[distance_metric_matrix==i]*=float(distribution[0][0])/distribution[0][i-1]
	#exp_dec=np.exp(-(distance_metric_matrix-1)/2.)
	return exp_dec

def transient_disposition_distribution(N,disp_distr_t,yb,kolm_smir,transition_flag):
	"""
	This function governs the transition from between two distributions 
	the functions similarity is derived via a kolmogorov smirnov test as a necessary criterium
	(https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test) 
	the target is set to 0.1, which is equivalent to the value one gets for random sampling from the given distribution

	The noise added is lognormal distributed (fat-tailed) to allow for long-range jumps and scaled with the deviation of the 
	actual distribution to the ideal curve. The greater the positive deviation, the greater the noise.

	"""
	kolm_smir_step=2
	improvement=True
	
	# Kolm Smir target 
	target=0.1
	# noise scaling coeff
	tr_noise_coeff=0.1

	def integrate_cdf(input_array):
		cdf=np.zeros(input_array.shape[0])
		for i in xrange(input_array.shape[0]):
			cdf[i]=integ.quad(lambda x: distribution_function(yb,x),0,input_array[i])[0]
		return cdf
	
	# set counters
	k=1
	n_inc=0
	
	# derive the ideal distribution for given yb
	x=np.linspace(0,1,101)
	target_dist=distribution_function(yb,x)
	target_dist=target_dist*float(N)/(sum(target_dist))
	# get deviation
	hist_values_sample,hist_bins=np.histogram(disp_distr_t,bins=100,range=(0,1))
	hist_values_sample=np.append(hist_values_sample,hist_values_sample[-1])
	hist_diff=hist_values_sample-target_dist
	# get the noise level for all N agents
	disp_distr_round=np.asarray(100*disp_distr_t,dtype='int')
	distr_dev=hist_diff[disp_distr_round]
	distr_dev=np.asarray(distr_dev,dtype='float')
	# scale and set positive
	distr_dev=distr_dev/50
	distr_dev+=np.abs(distr_dev.min())
	logn_std=1.

	if transition_flag== False:
		# incremental improvement for the non-transient case
		target=kolm_smir
	while kolm_smir_step >= target:
		# generate random lognormal dist and sign
		random_noise_increase=np.random.lognormal(0, logn_std, N)
		random_sign=np.random.randint(-1,1,N)
		random_sign[random_sign==0]=1
		# add the noise
		disp_distr_tp1=disp_distr_t+tr_noise_coeff*random_sign*random_noise_increase*distr_dev
		# check for boundaries
		disp_distr_tp1[disp_distr_tp1>1]=disp_distr_t[disp_distr_tp1>1]
		disp_distr_tp1[disp_distr_tp1<0]=disp_distr_t[disp_distr_tp1<0]
		# kolmogorov smirnov test
		kolm_smir_step=scipy.stats.kstest(disp_distr_tp1, integrate_cdf)[0]
		#print 'kolm_smir_step',kolm_smir_step	
		k+=1
		# in case of non-convergence, increase the standard deviation for the lognormal dist to allow for 
		if k>100:
			logn_std=logn_std*1.2
			print 'increased noise, logn_std',logn_std
			k=1
			n_inc+=1
			if n_inc>10:
				print 'DISTRIBUTION TRANSITION FAILED', 'kolm_smir_step',kolm_smir_step, 'yb',yb
				return disp_distr_t,kolm_smir_step,improvement
	return disp_distr_tp1,kolm_smir_step,improvement


def generate_initial_distance_sm(L):
	"""
	The initial proximity structure is generated. It's a small world NW following the 
	classical Watts & Strogatz approach with an initial ring that is successively rewired 
	with the rewiring probability L.p_rew. According to Watts & Strogatz 0.01 <p_rew < 0.1
	"""
	substr_adj_list=[]
	full_adjacency=np.zeros((L.N,L.N))
	proximity_small_world=np.zeros((L.N,L.N))
	proxim_nw=igraph.GraphBase.Lattice([L.N], nei=L.links_per_node, directed=False, mutual=True, circular=True)
	proxim_nw.rewire(int(L.p_rew*L.N*L.links_per_node))		
	small_world_distance_matrix=np.asarray(proxim_nw.shortest_paths())
	
	#derive the proximity matrix on the basis of the small-world network, random noise added to generate a continuum
	random_prox_noise=np.random.random((L.N,L.N))
	random_prox_noise=(random_prox_noise+random_prox_noise.transpose())*.5
	proximity_small_world=1-.1*(small_world_distance_matrix-1)-0.1*random_prox_noise

	# introduce a proximity offset,
	#if the distance falls below 5, there should be now considerable difference anymore
	proximity_small_world[np.where(proximity_small_world<=0.2)]=0.2
	
	# ensure that the diagonals are 0				
	for k in xrange(proximity_small_world.shape[0]):
		proximity_small_world[k,k]=0
		
	return proximity_small_world


def calc_cond_prob(smokers,nw_full,deg_sep_max,N):
	rcp=np.zeros(5)
	for i in xrange(deg_sep_max):
		deg_sep=i+1
		smoking_dep=[]
		for node in smokers:
			#print 'node',node
			distance_matrix=nw_full.path_lengths()	
			acquaintance_one=np.where(distance_matrix[node,:]==deg_sep)
			if acquaintance_one[0].size>0:
				smoking_dep.append(np.sum(nw_full.node_attribute('smoker')[acquaintance_one])/float(acquaintance_one[0].size)/(float(len(smokers))/N)-1)
		rcp[i]=np.mean(smoking_dep)
	return rcp

def derive_nw_chars(outdic,model_trans,kolm_smir_trans,L,i=0):
	if outdic.has_key('smokers')==False:
		l=1
		if L.transition_flag: l=L.n_iterations
		######################
		# initialize output dictionary
		######################
		outdic['no_interactions']=np.zeros(l)
		outdic['cogn_diss_coeff']=np.zeros(l)
		outdic['smokers']=[]
		outdic['no_of_smokers']=np.zeros(l)
		outdic['mean_degree']=np.zeros(l)
		outdic['clustering']=np.zeros(l)
		outdic['apl']=np.zeros(l)
		outdic['acquaintance_change']=np.zeros(l)
		outdic['centrality']={'smoker':np.zeros(l),'non_smoker':np.zeros(l)}
		if L.calc_full_centr_measures == True:
			outdic['betweenness']={'smoker':np.zeros(l),'non_smoker':np.zeros(l)}
			outdic['closeness']={'smoker':np.zeros(l),'non_smoker':np.zeros(l)}
			outdic['degree']={'smoker':np.zeros(l),'non_smoker':np.zeros(l)}
			outdic['ind_clustering']={'smoker':np.zeros(l),'non_smoker':np.zeros(l)}
		outdic['no_disconnected_nodes']=np.zeros(l)
		outdic['kolm_smir']=np.zeros(l)
		outdic['conditional_prob']=np.zeros((l,5))
		outdic['disp_distr']=np.zeros((l,L.N))
		outdic['max_sm_cl_size']=np.zeros(l)
		outdic['ave_sm_cl_size']=np.zeros(l)

	########################################################################################################
	# write output model chars
	outdic['no_interactions'][i]=model_trans._no_interactions
	# print 'no_interactions',outdic['no_interactions'][i]
	outdic['acquaintance_change'][i]=model_trans._new_edges
	smokers=np.where(model_trans.get_char_distribution()[1,:]==1)
	non_smokers=np.where(model_trans.get_char_distribution()[1,:]==0)
	outdic['smokers'].append(smokers)
	outdic['no_of_smokers'][i]=len(outdic['smokers'][-1][0])
	outdic['cogn_diss_coeff'][i]=np.mean(np.abs(model_trans.get_char_distribution()[0,:]-model_trans.get_char_distribution()[1,:]))
	outdic['disp_distr'][i,:]=model_trans.get_char_distribution()[0,:]

	#######
	# derive max smoker cluster size
	ts_adj=model_trans.get_acquaintance_network().adjacency
	sm_mat=np.zeros(1000)
	sm_mat[smokers]=1
	smok=np.asarray(sm_mat,dtype='bool')
	try:
		nwsm=Network(adjacency=ts_adj[smok,:][:,smok])
		get_cl=nwsm.graph.clusters(mode='WEAK')
		sm_cl_len=[ len(get_cl[rft]) for rft in xrange(len(get_cl))]
		sm_cl_len=np.asarray(sm_cl_len)
		outdic['max_sm_cl_size'][i]=np.max(sm_cl_len)
		# print sm_cl_len[sm_cl_len>1]
		outdic['ave_sm_cl_size'][i]=np.mean(sm_cl_len[sm_cl_len>1])
	except: print 'conditional prob calc failed'

	#write output network chars
	outdic['mean_degree'][i]=model_trans.get_acquaintance_network().degree().mean()
	outdic['clustering'][i]=model_trans.get_acquaintance_network().global_clustering()
	outdic['kolm_smir'][i]=kolm_smir_trans
	outdic['apl'][i]=model_trans.get_acquaintance_network().average_path_length()
	outdic['conditional_prob'][i,:]=calc_cond_prob(smokers[0],model_trans.get_acquaintance_network(),L.cond_prob_degree,L.N)
	#print 'apl', outdic['apl'][i]
	outdic['centrality']['smoker'][i]=np.mean(np.asarray(model_trans.get_acquaintance_network().graph.evcent(scale=False))[smokers])
	outdic['centrality']['non_smoker'][i]=np.mean(np.asarray(model_trans.get_acquaintance_network().graph.evcent(scale=False))[non_smokers])

	if L.calc_full_centr_measures == True:
		outdic['betweenness']['smoker'][i]=np.mean(model_trans.get_acquaintance_network().betweenness()[smokers])
		outdic['betweenness']['non_smoker'][i]=np.mean(model_trans.get_acquaintance_network().betweenness()[non_smokers])
		outdic['closeness']['smoker'][i]=np.mean(model_trans.get_acquaintance_network().closeness()[smokers])
		outdic['closeness']['non_smoker'][i]=np.mean(model_trans.get_acquaintance_network().closeness()[non_smokers])
		degree_matrix=model_trans.get_acquaintance_network().degree()
		outdic['degree']['smoker'][i]=np.mean(degree_matrix[smokers])
		outdic['degree']['non_smoker'][i]=np.mean(degree_matrix[non_smokers])
		outdic['ind_clustering']['smoker'][i]=np.mean(model_trans.get_acquaintance_network().local_clustering()[smokers])
		outdic['ind_clustering']['non_smoker'][i]=np.mean(model_trans.get_acquaintance_network().local_clustering()[non_smokers])
		disc_node=(np.where(np.isinf(model_trans.get_acquaintance_network().path_lengths()[1,:])==True)[0])
		outdic['no_disconnected_nodes'][i]=disc_node.shape[0]

	# print i, 'i'
	# print 'apl',outdic['apl'][i]
	# print 'no_of_smokers', outdic['no_of_smokers'][i]
		
	
	return outdic

####################################################################################################################################
# SINGLE RUN
####################################################################################################################################
def generate_eq(L):
	print '############## INITIALIZE MODEL'
	no_components=2
	while no_components !=1:

		char_distribution_initial= rejection_sampling(L.N, distribution_function,3.)

		#  Degree preference distribution
		degree_preference = np.random.normal(L.mean_degree_pref, L.std_degree_pref, L.N).astype("int8")
		
		# Get the underlying network structure
		proximity_structure=generate_initial_distance_sm(L)
				
		###################
		# Generate equilibrium acqaintance networks
		###################
		
		#  Create random acquaintance network (Erdös-Renyi graph)
		acquaintance_network = Network.ErdosRenyi(n_nodes=L.N, link_probability=L.mean_degree_pref / float(L.N - 1))
		
		#  Set silence level
		acquaintance_network.silence_level = 3
		
		char_feedback=True
		dynamic=True
		model_initial=DynamicSocietyModel (proximity_structure,char_distribution_initial, acquaintance_network,
	                 interaction_likelihood_function,L,char_feedback,dynamic,degree_preference)
		model_initial.get_eq_network()
		model_initial.get_acquaintance_network().set_node_attribute('smoker',model_initial.get_char_distribution()[1,:])

		print '############## RUN EQ GENERATION'
		for i in xrange(L.n_initial_eq_it):
			model_initial.iterate(1)
			model_initial.get_acquaintance_network().set_node_attribute('smoker',model_initial.get_char_distribution()[1,:])
		print '--------------- EQ Network generated --------------- '

		clusters=model_initial.get_acquaintance_network().graph.clusters(mode='WEAK')
		# len(model_hysteresis.get_acquaintance_network().graph.clusters())
		no_components= len(clusters[0])+len(clusters)-L.N
		print '################# Number of components initial network #################'
		print no_components,len(clusters)

		# check for feasibility of the transient disp distr
		# if kolm_smir_trans initial> 0.1 than the procedure fails
		kolm_smir_trans=0.
		char_dist=model_initial.get_char_distribution()[0,:]
		char_dist,kolm_smir_trans,improvement_random_process=transient_disposition_distribution(L.N,char_dist,3.,kolm_smir_trans,L.transition_flag)
		if kolm_smir_trans>0.1: 
			print 'initial transition test failed'
			no_components+=10**6
	return model_initial

def transition(char_feedback,dynamic,model_initial,L,out):
	model_trans = DynamicSocietyModel(model_initial._proximity_structure,model_initial.get_char_distribution(), model_initial.get_acquaintance_network(),
                  interaction_likelihood_function,L,char_feedback,dynamic,model_initial._degree_preference)	
	######################
	# run_script
	######################
	if 	char_feedback==True:
		if dynamic== True:
			print '############## RUN FULLY COUPLED'
			nw_snapshots_dic={}
		else: print '############## RUN SOCIAL INFLUENCE ONLY'
	else: 
		if L.dyn_meanfield:
			print '############## RUN MEAN FIELD FORCING'
		else:
			print '############## RUN DYNAMIC ONLY'
	
	######################
	# derive the final distribution
	######################
		
	kolm_smir_trans=1
	improvement_random_process=True
	print 'transient part'
	model_trans.get_acquaintance_network().set_node_attribute('smoker',model_trans.get_char_distribution()[1,:])
	output_dict={}
	nw_snapshots_dic={}
	for i in xrange(L.n_iterations):
		
		###############################################################	
		#	TRANSIENT CHANGE OF THE SMOKING DISPOSITION
		###############################################################
		if L.transition_flag==True:
			if i<L.n_transition:
				yb=3.-(3.-L.yb_final)*(float(i+1)/L.n_transition)
				# print 'yb',yb,i
				char_dist=model_trans.get_char_distribution()[0,:]
				char_dist_step,kolm_smir_trans,improvement_random_process=transient_disposition_distribution(L.N,char_dist,yb,kolm_smir_trans,L.transition_flag)
				model_trans.set_char_distribution(char_dist_step,0)
		
				if char_feedback==False and L.dyn_meanfield==False:
					model_trans.set_char_distribution(np.round(model_trans.get_char_distribution()[0,:],0),1)
					#print np.round(model_trans.get_char_distribution()[0,:10],0)
			
			##############################################################
			model_trans.iterate(1)
			model_trans.get_acquaintance_network().set_node_attribute('smoker',model_trans.get_char_distribution()[1,:].copy())
			#print 'no of smokers',sum(model_trans.get_char_distribution()[1,:])
			output_dict=derive_nw_chars(output_dict,model_trans,kolm_smir_trans,L,i)

		###############################################################	
		#	ABRUPT CHANGE OF THE SMOKING DISPOSITION
		###############################################################
		else:
			char_dist=model_trans.get_char_distribution()[0,:]
			if kolm_smir_trans>0.1:
				char_dist_step,kolm_smir_trans,improvement_random_process=transient_disposition_distribution(L.N,char_dist,L.yb_final,kolm_smir_trans,L.transition_flag)
				model_trans.set_char_distribution(char_dist_step,0)
			if char_feedback==False and L.dyn_meanfield==False:
				model_trans.set_char_distribution(np.round(model_trans.get_char_distribution()[0,:],0),1)
			print 'no of smokers',sum(model_trans.get_char_distribution()[1,:])
			print 'kolm_smir_trans',kolm_smir_trans
			##############################################################
			model_trans.iterate(1)
			model_trans.get_acquaintance_network().set_node_attribute('smoker',model_trans.get_char_distribution()[1,:].copy())
		#######################################################################################################
		# SAVE NETWORK INSTANCES
		#######################################################################################################

		if (char_feedback&dynamic==True) and (i==L.n_iterations-1 or i in np.arange(0,L.n_iterations,L.nw_save_steps,dtype='int')):
			nw_snapshots_dic[i]=(model_trans.get_acquaintance_network().adjacency)

	if (char_feedback&dynamic==True):
		if L.transition_flag: output = open(L.output_path+'/nw_snaps_trans_smok_%s.pkl'%(out), 'wb')
		else: output = open(L.output_path+'/nw_snaps_eq_yb_%s_%s.pkl'%(L.yb_final,out), 'wb')
		pickle.dump(nw_snapshots_dic, output)
		output.close()

	if L.transition_flag==True:
		return output_dict
	else: 
		return derive_nw_chars(output_dict,model_trans,kolm_smir_trans,L)



def do_one(out,L):
	print 'out', out
	####################################################################################################################################
	#
	#  Main script
	#
	####################################################################################################################################
	
	model_initial= generate_eq(L)
	out_dic_2x2={}
	
	if L.coupling_instances['full']:
		out_dic_2x2['full']=transition(True,True,model_initial,L,out)
	if L.coupling_instances['infl']:
		out_dic_2x2['infl']=transition(True,False,model_initial,L,out)
	if L.coupling_instances['dyn']:
		# Flag setting the mean field for dynamical forcing, historically coupled to the dyn case
		L.dyn_meanfield=False
		out_dic_2x2['dyn']=transition(False,True,model_initial,L,out)
	if L.coupling_instances['mean_field']:
		# Flag setting the mean field for dynamical forcing, historically coupled to the dyn case
		L.dyn_meanfield=True
		out_dic_2x2['mean_field']=transition(False,True,model_initial,L,out)
		
	if L.transition_flag: output = open(L.output_path+'/trans_smok_%s.pkl'%(out), 'wb')
	else: output = open(L.output_path+'/eq_yb_%s_%s.pkl'%(L.yb_final,out), 'wb')
	pickle.dump(out_dic_2x2, output)
	output.close()

	
