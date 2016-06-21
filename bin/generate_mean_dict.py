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
path =''#'/scratch/01/carls/social/emulation/ens_rc_p25_ens_100/'

run_reps=glob.glob('trans_smok*')
ensemble_size=len(run_reps)
summary_dic={}

for run in run_reps:
	pkl_file=open(run,'rb')
	read_in=pickle.load(pkl_file)
	pkl_file.close()
	print run
	for active_key in read_in.iterkeys():

		for active_var in read_in[active_key].iterkeys():
		# for active_var in read_out_list:
			#print active_var, active_ke
			if summary_dic.has_key(active_var) == False:
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

	#centr_connected.append((folder,np.asarray(no_connected_list).mean(),np.asarray(no_connected_list).std()))

#############################
# SAVE OUTPUT
#############################

output = open('full_trans_smok.pkl', 'wb')
pickle.dump(summary_dic, output)
output.close()
