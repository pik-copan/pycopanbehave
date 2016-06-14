#!/usr/bin/python
# -*- coding: utf-8 -*-

# Add repository path of the model core
import sys
import os
import glob
sys.path.append('/Users/carls/Documents/PIK/git_repos/pynamic-society/bin')
sys.path.append('/home/carls/git_repos/pynamic-society/bin')

from ens_smoking_transition import *

# output path

#output_path='/scratch/01/carls/social/mod_smok_mobility/smok_mob_1'

output_path='./'#/Users/carls/Downloads/social/base_proxim_02_rand_05_sm_mob_2'

#check if output path exists
if os.path.exists(output_path)==False:
	sys.exit("FATAL: output directory does not exist")

#######################################
# Namelist as a Bunch dictionary
L=Bunch()

L.output_path=output_path

#  Ensemble size (number of realizations of model time evolution)
L.n_ensemble = 96

# Flag, if a transition from the initial to the final distribution should take place. If false, the final distribution will be applied directly
L.transition_flag=True

# Flag, which instances should be run (full,infl,dyn)
L.coupling_instances={'full':True,'infl':True,'dyn':True}

# Flag, if the dyn only case should be derived via a mean-field
L.dyn_meanfield=True

# sets the final value for the disposition function (yb_final==0.0 is full transition)
L.yb_final=0.5

#  Number of nodes
L.N = 1000

#  Mean degree preference
L.mean_degree_pref = 10 # as it is stated in cf08
#  Degree preference standard deviation
L.std_degree_pref =3

# Rewiring probability of the underlying proximity 
L.p_rew=0.03
L.links_per_node=2

#probability of interaction
L.p_ai=.8

# offset that sets a basic interaction probability with all agents
L.interaction_offset=0.03

# sets the mobility with regard to the underlying proximity structure in positions
smoking_mobility=2
smoking_weight=.1*smoking_mobility
L.char_weight=(0,smoking_weight,(1-smoking_weight))

#  Number of hysteresis iterations
L.n_transition = 1000
L.n_initial_eq_it=200
L.n_trans_eq_post=0
L.n_iterations=L.n_transition+L.n_trans_eq_post

"""
Parameters only relevant for the transition_flag==True case 
"""

# Flag, if betweeneess and closeness should be recorded calculated transient online (slow)
L.calc_full_centr_measures=True

# Sets the level of degree up to which conditional probability should be derived (max 5)
# warning, the larger, the slower
L.cond_prob_degree=5

# number of snapshots of the full nw
nw_snapshots=10
L.nw_save_steps=int(L.n_iterations/nw_snapshots)

"""
Parameters only relevant for the transition_flag==False case 
"""
discrete_stops_over_interval=20

yb_initial=3.0

####################################################################################################################################
#
#  MPI script
#
####################################################################################################################################


def master():
	out=0
	if L.transition_flag:
		for i in xrange(L.n_ensemble):
			mpi.submit_call("do_one",(i,L),id=i)
			out+=1
	else:
		paramrange=np.linspace(yb_initial,L.yb_final,discrete_stops_over_interval)
		# for k in xrange(discrete_stops_over_interval):	
		# 	for i in xrange(L.n_ensemble):
		L.yb_final=3.0#paramrange[k]
		mpi.submit_call("do_one",(out,L),id=out)
		out+=1	

mpi.run()


# no_runs=len(glob.glob('nw_snaps_trans_smok*pkl'))
# for run_no in xrange(no_runs):
# 	#run_no=str(0)
# 	nwdi=np.load('nw_snaps_trans_smok_'+str(run_no)+'.pkl')
# 	prop_dict=np.load('trans_smok_'+str(run_no)+'.pkl')
# 	for ts in (0,900):#nwdi.iterkeys():
# 		# ts=0
# 		sm_mat=np.zeros(1000)
# 		smokers=prop_dict['full']['smokers'][ts]
# 		sm_mat[smokers]=1
# 		nw=Network(adjacency=nwdi[ts])
# 		nw.set_node_attribute('smoker',sm_mat)
# 		nw.save('nw_snapshot_run_'+str(run_no)+'_timestep_'+str(ts)+'.graphml',format='graphml')


