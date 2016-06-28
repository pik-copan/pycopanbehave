#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
pycopanbehave -- An adaptive network mode of behaviour selection in Python

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

from ens_smoking_transition import *

sys.path.append('../bin')


#
#  Initializations
#

output_path = './ens_members/'

# check if output path exists
if not os.path.exists(output_path):
    os.mkdir(output_path)

#######################################
# Namelist as a Bunch dictionary
L = Bunch()

L.output_path = output_path

#  Ensemble size (number of realizations of model time evolution)
L.n_ensemble = 1000

# Flag, if a transition from the initial to the final distribution should take
# place. If false, the final distribution will be applied directly
L.transition_flag = True

# Flag, which instances should be run (full,infl,dyn,meanfield)
L.coupling_instances = {'full': True, 'infl': True, 'dyn': True,
                        'mean_field': True}

# Sets the final value for the disposition function (yb_final==0.0 is full
# transition)
L.yb_final = 0.5

#  Number of nodes
L.N = 100

#  Mean degree preference
L.mean_degree_pref = 10  # as it is stated in cf08
#  Degree preference standard deviation
L.std_degree_pref = 3

# Rewiring probability of the underlying proximity
L.p_rew = 0.03
L.links_per_node = 2

# Probability of interaction
L.p_ai = .8

# Offset that sets a basic interaction probability with all agents
L.interaction_offset = 0.03

# Sets the mobility with regard to the underlying proximity structure in
# positions
smoking_mobility = 2
smoking_weight = .1 * smoking_mobility
L.char_weight = (0, smoking_weight, (1 - smoking_weight))

# Smoking behaviour switching probability scaling factor  scales the switching
# probability pi(t) of the smoking
# behaviour.
# C controls the amplitude of equilibrium stochastic noise of the smoking
# behaviour that is introduced by the Ising-like implementation.
L.C = 0.1

#  Number of hysteresis iterations
L.n_transition = 1000
L.n_initial_eq_it = 200
L.n_trans_eq_post = 0
L.n_iterations = L.n_transition + L.n_trans_eq_post

#
#  Parameters only relevant for the transition_flag==True case
#

# Flag, if betweeneess and closeness should be recorded calculated transient
# online (slow)
L.calc_full_centr_measures = True

# Sets the level of degree up to which conditional probability should be
# derived (max 5)
# warning, the larger, the slower
L.cond_prob_degree = 5

# number of snapshots of the full nw
nw_snapshots = 10
L.nw_save_steps = int(L.n_iterations / nw_snapshots)

#
#  Parameters only relevant for the transition_flag==False case
#

discrete_stops_over_interval = 20

yb_initial = 3.0

###############################################################################
#
#  MPI script
#
###############################################################################


def master():
    out = 0

    if L.transition_flag:
        for i in xrange(L.n_ensemble):
            mpi.submit_call("do_one", (i, L), id=i)
            out += 1
    else:
        paramrange = np.linspace(yb_initial, L.yb_final
                                 discrete_stops_over_interval)
        L.yb_final = 3.0  # paramrange[k]
        mpi.submit_call("do_one", (out, L), id=out)
        out += 1

mpi.run()
