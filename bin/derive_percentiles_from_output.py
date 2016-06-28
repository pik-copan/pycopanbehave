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

SUPPORT FUNCTION TO DERIVE PERCENTILES FROM INDIVDUAL MODEL OUTPUT
"""

#
#  Imports
#

import numpy as np
import glob
import pickle
import scipy.stats

range = getattr(__builtins__, 'xrange', __builtins__.range)

#
#  Initializations
#

# SET INPUT PATH OF ENSEMBLE SIMULATIONS
path = 'ens_members/'


run_reps = glob.glob(path + 'trans_smok*')
ensemble_size = len(run_reps)
ens_dic = {}
full_dict = {}
dict_disconnected = {}

pkl_file = open(run_reps[0], 'rb')
read_in = pickle.load(pkl_file)
pkl_file.close()
fin_dic = {}


for active_var in read_in[read_in.keys()[0]].keys():
    fin_dic[active_var] = {}

summary_dic = {}
no_connected_list = []
folder = 'ens'
folder_key = str(folder)[-2:]
dict_disconnected[folder] = {}
for run in range(ensemble_size):
    pkl_file = open(path + 'trans_smok_%s.pkl' % run, 'rb')
    read_in = pickle.load(pkl_file)
    print(run)
    for active_key in read_in.keys():

        for active_var in read_in[active_key].keys():
            if active_var == 'centrality_smokers_connected':
                if active_key == 'full':
                    rel = read_in['full']['centrality']['smoker'][:10].mean()
                    no_connected_list.append(
                        read_in['full']['centrality_smokers_connected']/rel)

            elif active_var == 'disconnected_nodes':
                if active_key == 'full':
                    dict_disconnected[folder][run] = \
                        read_in[active_key][active_var]

            elif active_var not in summary_dic:
                summary_dic[active_var] = {}

            elif active_var == 'smokers':
                if active_key not in summary_dic[active_var]:
                    summary_dic[active_var][active_key] = []
                summary_dic[active_var][active_key].append(
                    read_in[active_key][active_var])

            elif np.asarray(read_in[active_key][active_var]).dtype == 'object':
                for sm_key in read_in[active_key][active_var].keys():
                    if sm_key not in summary_dic[active_var]:
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

    pkl_file.close()

#############################
# DERIVE ENSEMBLE PERCENTILES
#############################

print('DERIVE ENSEMBLE PERCENTILES')

ave_dic = {}
for active_var in summary_dic.keys():
    if active_var not in ave_dic:
        ave_dic[active_var] = {}

    print(active_var, active_key)
    if active_var == 'no_disconnected_nodes':
        for run_key in summary_dic[active_var].keys():
            if run_key not in ave_dic[active_var]:
                ave_dic[active_var][run_key] = {}
            dataset = np.asarray(summary_dic[active_var][run_key])
            dataset[dataset > 100] = 1
            ave_dic[active_var][run_key]['mean'] = dataset.mean(axis=0)
            ave_dic[active_var][run_key]['lower'] = dataset.min(axis=0)
            ave_dic[active_var][run_key]['upper'] = dataset.max(axis=0)

    elif 'smoker' in summary_dic[active_var]:
        print(active_var)
        for active_key in summary_dic[active_var].keys():
            if active_key not in ave_dic[active_var]:
                ave_dic[active_var][active_key] = {}

            for run_key in summary_dic[active_var][active_key].keys():
                if run_key not in ave_dic[active_var][active_key]:
                    ave_dic[active_var][active_key][run_key] = {}
                    dataset = np.asarray(
                        summary_dic[active_var][active_key][run_key])
                    ave_dic[active_var][active_key][run_key]['mean'] \
                        = scipy.stats.scoreatpercentile(
                            dataset / np.mean(dataset[:, 1:25], axis=1)
                                        .reshape(dataset.shape[0], 1), 50)
                    ave_dic[active_var][active_key][run_key]['lower'] \
                        = scipy.stats.scoreatpercentile(
                            dataset /
                            np.mean(dataset[:, 1:25], axis=1).reshape(
                                dataset.shape[0], 1), 25)
                    ave_dic[active_var][active_key][run_key]['upper'] \
                        = scipy.stats.scoreatpercentile(
                            dataset / np.mean(dataset[:, 1:25], axis=1)
                                        .reshape(dataset.shape[0], 1), 75)
    elif np.asarray(summary_dic[active_var]['full']).dtype == object:
        print('smokers yes')

        for run_key in summary_dic[active_var].keys():
            if run_key not in ave_dic[active_var]:
                ave_dic[active_var][run_key] = {}
            dataset = np.zeros((len(summary_dic[active_var][run_key]),
                                len(summary_dic[active_var][run_key][0])))
            print("dataset.shape", dataset.shape)
            for i in range(len(summary_dic[active_var][run_key])):
                for k in range(len(summary_dic[active_var][run_key][i])):
                    x = summary_dic[active_var][run_key][i][k][0].shape[0]
                    dataset[i, k] = x

            ave_dic[active_var][run_key]['mean'] = dataset.mean(axis=0)
            ave_dic[active_var][run_key]['lower'] = dataset.min(axis=0)
            ave_dic[active_var][run_key]['upper'] = dataset.max(axis=0)

    else:
        for run_key in summary_dic[active_var].keys():
            if run_key not in ave_dic[active_var]:
                ave_dic[active_var][run_key] = {}
            dataset = np.asarray(summary_dic[active_var][run_key])
            ave_dic[active_var][run_key]['mean'] \
                = scipy.stats.scoreatpercentile(dataset, 50)
            ave_dic[active_var][run_key]['lower'] \
                = scipy.stats.scoreatpercentile(dataset, 25)
            ave_dic[active_var][run_key]['upper'] \
                = scipy.stats.scoreatpercentile(dataset, 75)

#############################
# DERIVE FINAL VALUES
#############################

# print(runsumdict)
run_sum_dict = summary_dic  # XXX bug? this will be no copy
for active_var in run_sum_dict.keys():
    if 'smoker' in run_sum_dict[active_var]:
        print(active_var)
        for active_key in run_sum_dict[active_var].keys():
            dataset = np.asarray(run_sum_dict[active_var]['smoker']['full'])
            fin_dic[active_var][folder_key] = np.mean(
                dataset, axis=1) / np.mean(dataset[:, 1:25], axis=1)
    elif np.asarray(run_sum_dict[active_var]['full']).dtype == object:
        print('smokers yes')
        for run_key in run_sum_dict[active_var].keys():
            dataset = np.zeros((len(run_sum_dict[active_var][run_key]),
                                len(run_sum_dict[active_var][run_key][0])))
            print("dataset.shape", dataset.shape)
            for i in range(len(run_sum_dict[active_var][run_key])):
                for k in range(len(run_sum_dict[active_var][run_key][i])):
                    dataset[i, k] = \
                        run_sum_dict[active_var][run_key][i][k][0].shape[0]
            fin_dic[active_var][folder_key] = np.mean(dataset, axis=1)

    else:
        for run_key in run_sum_dict[active_var].keys():
            if run_key not in fin_dic[active_var]:
                fin_dic[active_var][run_key] = {}
            dataset = np.asarray(run_sum_dict[active_var][run_key])
            fin_dic[active_var][folder_key] = np.mean(dataset, axis=1)

#############################
# SAVE OUTPUT
#############################

with open('median_trans_smok.pkl', 'wb') as fid:
    pickle.dump(ave_dic, fid)
ens_dic = ave_dic
full_dict = summary_dic

with open('full_trans_smok.pkl', 'wb') as fid:
    pickle.dump(full_dict, fid)

with open('ens_trans_smok.pkl', 'wb') as fid:
    pickle.dump(ens_dic, fid)

with open('final_reduction.pkl', 'wb') as fid:
    pickle.dump(fin_dic, fid)
