#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Smoking Transition -- A use case of the COPAN:BEHAVE model to emulate
changes in smoking behaviour

Copyright (C) 2011--2016 Potsdam Institute for Climate Impact Research
Authors: Jonathan F. Donges <donges@pik-potsdam.de>,
         Carl-Friedrich Schleussner <schleussner@pik-potsdam.de>,
         Denis Engemann <denis.engemann@gmail.com>
URL:     <http://www.pik-potsdam.de/copan/software>

RUN FILE SMOKER TRANSITION
"""

#
#  Imports
#

# Add repository path of the model core
import sys
import os
import glob
sys.path.append('../bin')

from ens_smoking_transition import *

#
#  Initializations
#

output_path = './ens_members/'

# Check if output path exists
if os.path.exists(output_path) == False:
	os.mkdir(output_path)

#######################################
# Namelist as a Bunch dictionary
L = Bunch()

#
#  Parameter settings
#

L.output_path = output_path

# Ensemble size (number of realizations of model time evolution)
L.n_ensemble = 3

# Flag, if a transition from the initial to the final distribution of
# smoking disposition should take place. If false, the final distribution
# will be applied directly
L.transition_flag = True

# Flag, which instances should be run (full, infl, dyn, meanfield)
L.coupling_instances = {'full':True, 'infl':True, 'dyn':True,
						'mean_field':True}

# Placeholder flag for mean field dynamics
L.dyn_meanfield = None

# Sets the final value for the smoking disposition function
# (yb_final = 0.0 is full transition)
L.yb_final = 1.5

#  Number of nodes
L.N = 500

# Mean degree preference (value as stated in Christakis and Fowler, The
# New England Journal of Medicine 358, 2008)
L.mean_degree_pref = 10
# Degree preference standard deviation
L.std_degree_pref = 3

# Rewiring probability of the underlying proximity
L.p_rew = 0.03
L.links_per_node = 2

# Probability of interaction
L.p_ai = .8

# Offset that sets a basic interaction probability with all agents
L.interaction_offset = 0.03

# Sets the mobility with regard to the underlying proximity structure
# in positions
smoking_mobility = 2
smoking_weight = .1 * smoking_mobility
L.char_weight = (0, smoking_weight, (1 - smoking_weight))

#  Number of hysteresis iterations
L.n_transition = 100
L.n_initial_eq_it = 10
L.n_trans_eq_post = 0
L.n_iterations = L.n_transition + L.n_trans_eq_post

#
# Parameters only relevant for the transition_flag==True case
#

# Flag, if betweeneess and closeness should be recorded calculated transient online (slower)
L.calc_full_centr_measures = True

# Sets the level of degree up to which conditional probability should be derived (max 5, the larger, the slower)
L.cond_prob_degree = 5

# number of snapshots of the full nw
nw_snapshots = 2
L.nw_save_steps = int(L.n_iterations / nw_snapshots)

#
# Parameters only relevant for the transition_flag==False case
#

discrete_stops_over_interval = 20

yb_initial = 3.0

####################################################################################################################################
#
#  MPI script
#
####################################################################################################################################


def master():
	out = 0
	if L.transition_flag:
		for i in xrange(L.n_ensemble):
			mpi.submit_call("do_one",(i,L),id=i)
			out += 1
	else:
		paramrange=np.linspace(yb_initial,L.yb_final,discrete_stops_over_interval)
		# for k in xrange(discrete_stops_over_interval):
		# 	for i in xrange(L.n_ensemble):
		L.yb_final = 3.0 #paramrange[k]
		mpi.submit_call("do_one",(out,L),id=out)
		out += 1

mpi.run()