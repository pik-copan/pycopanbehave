#!/usr/bin/python
# -*- coding: utf-8 -*-

#
#  Imports
#

#  Import numpy for numerics
import numpy as np


# Add repository path of the model core
import sys

sys.path.append('/Users/carls/Documents/PIK/ModelingSocialStability/Code/dynamic_society_model/src/model_core')
sys.path.append('/Users/carls/Documents/PIK/ModelingSocialStability/Code/src/model_core')
sys.path.append('/scratch/01/carls/social/emulation/model_core/')
sys.path.append('../../../')
sys.path.append('../../')
#  Import Network class from pygeonetwork
# from pygeonetwork import Network

#  Import DynamicSocietyModel class from pydynamicsociety
# from pysoc import DynamicSocietyModel

# from pygeonetwork import mpi

import pickle

import scipy.stats
import scipy.integrate as integ



#######################################
#  Number of nodes
N = 1000

#  Number of hysteresis iterations
n_transition = 1000


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


char_dist=np.zeros((N,n_transition))
char_dist[:,0]= rejection_sampling(N, distribution_function,3)[0,:]
# char_dist=model_trans._char_distribution[0,:]
kolm_smir_trans=2
for i in xrange(1,n_transition):
	yb_final=3.0
	yb=3.-yb_final*(float(i+1)/n_transition)
	print 'yb',yb,i						
	char_dist[:,i],kolm_smir_trans,improvement_random_process=transient_disposition_distribution(N,char_dist[:,i-1],yb,kolm_smir_trans,True)


#kolm_smir[i]=kolm_smir_trans

yb=1.83
tr_noise_coeff=0.1*1000/n_transition
disp_distr_t=rejection_sampling(N, distribution_function,3)[0,:]
sample_dist=rejection_sampling(N, distribution_function,yb)[0,:]
hist_values_sample,hist_bins=np.histogram(disp_distr_t,bins=100,range=(0,1))
hist_values_ref,hist_bins=np.histogram(rejection_sampling(N, distribution_function,yb)[0,:],bins=100,range=(0,1))

hist_diff=hist_values_ref-hist_values_sample
hist_diff=np.append(hist_diff,hist_diff[-1])
disp_distr_round=np.asarray(100*disp_distr_t,dtype='int')

distr_dev=hist_diff[disp_distr_round]
distr_dev=np.asarray(distr_dev,dtype='float')
distr_dev=distr_dev/50#np.abs(distr_dev).max()

random_noise_increase=np.random.rand(N) 

disp_distr_tp1=tr_noise_coeff*random_noise_increase*distr_dev


point_stability=distribution_function(yb,disp_distr_t)
point_stability=1-point_stability/point_stability.max()
disp_distr_tp1=tr_noise_coeff*(np.random.randn(N))*point_stability

plot(disp_distr_t,disp_distr_tp1,'x')
