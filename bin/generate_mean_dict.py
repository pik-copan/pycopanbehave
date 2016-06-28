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

SUPPORT FUNCTIONS
"""

#
#  Imports
#

import sys
import glob
import os
import pickle

import numpy.ma as ma
import numpy as np
import scipy.stats

#
#  Initializations
#

# SET INPUT PATH OF ENSEMBLE SIMULATIONS
path = ''  # '/scratch/01/carls/social/emulation/ens_rc_p25_ens_100/'

run_reps = glob.glob('trans_smok*')
ensemble_size = len(run_reps)
summary_dic = {}

for run in run_reps:
    with open(run, 'rb') as fid:
        read_in = pickle.load(fid)
    print(run)
    for active_key in read_in.keys():
        # for active_var in read_out_list:
        for active_var in read_in[active_key].keys():
            # print(active_var, active_ke)
            if active_var not in summary_dic:
                summary_dic[active_var] = {}

            elif active_var == 'smokers':
                # print('smokers unpacked')
                if active_key not in summary_dic[active_var]:
                    summary_dic[active_var][active_key] = []
                summary_dic[active_var][active_key].append(
                    read_in[active_key][active_var])

            elif np.asarray(read_in[active_key][active_var]).dtype == 'object':
                for sm_key in read_in[active_key][active_var].keys():
                    if if sm_key not insummary_dic[active_var]:
                        summary_dic[active_var][sm_key] = {}

                    if active_key not in summary_dic[active_var][sm_key]:
                        summary_dic[active_var][sm_key][active_key] = []

                    summary_dic[active_var][sm_key][active_key].append(
                        read_in[active_key][active_var][sm_key])
            else:
                if active_key not in summary_dic[active_var]:
                    summary_dic[active_var][active_key] = []

                summary_dic[active_var][active_key].append(
                    read_in[active_key][active_var])

    # centr_connected.append((folder,np.asarray(no_connected_list).mean(),
    # np.asarray(no_connected_list).std()))

#############################
# SAVE OUTPUT
#############################

with open('full_trans_smok.pkl', 'wb') as fid:
    pickle.dump(summary_dic, fid)
