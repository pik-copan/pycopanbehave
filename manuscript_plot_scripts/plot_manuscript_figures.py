#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
pycopanbehave -- An adaptive network mode of behaviour selection in Python

Copyright (C) 2011--2016 Potsdam Institute for Climate Impact Research
Authors: Jonathan F. Donges <donges@pik-potsdam.de>,
         Carl-Friedrich Schleussner <schleussner@pik-potsdam.de>,
         Denis Engemann <denis.engemann@gmail.com>
URL:     <http://www.pik-potsdam.de/copan/software>

PLOT SCRIPTS
"""

#
#  Imports
#

import os.path as op
import os
import pickle

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#
#  Initializations
#

input_file = 'full_trans_smok.pkl'
# input_file='trans_smok_output_full.pkl'
with open(input_file) as w:
    data = pickle.load(w)

if not os.path.isdir('figures'):
    os.mkdir('figures')

to_plot = [1, 2, 3]
set_titles = False
smokers = data['no_of_smokers']
decim = 1
ci = 66
n_smokers = 1000.
start_time = 0.
n_xticks = [0, 200, 400, 600, 790]
ma_window_size = 30
keys = ['full', 'dyn', 'infl', 'mean_field']
mapping = dict(zip(keys, ['coupled', 'network', 'interaction', 'mean-field']))
colors_context = "#55247A", "#3368A5", "#D33C3E", '#37442A'
color_smoker = '#A44791'
color_nonsmoker = '#BECD00'
sns.set(style='ticks', context='poster', font_scale=1.4)
my_gray = '#6e6a6a'

plt.rcParams.update({k: my_gray for k in (
    'axes.labelcolor',
    'axes.edgecolor',
    'text.color',
    'xtick.color',
    'ytick.color'
)})

plt.rcParams.update({k: 2.2 for k in (
    'axes.linewidth',
    'ytick.major.width',
    'xtick.major.width'
)})

#
#  DEFINE PLOTTING FUNCTIONS
#


def moving_average(values, window, axis):
    """
    Add docstring!
    """
    weights = np.repeat(1.0, window) / window
    sma = np.apply_along_axis(
        lambda m: np.convolve(m, weights, mode='valid'), arr=values, axis=axis)
    return sma


def preprocess_array(X):
    """
    Add docstring!
    """
    return moving_average(np.transpose(X, (1, 2, 0)), ma_window_size, 1)


def rescale(x, axis, method='percent'):
    """
    Add docstring!
    """
    if method == 'divide':
        x = x / x[:, 0:2, :]
    elif method == 'percent':
        x = (x - x[:, 0:1, :]) / x[:, 0:1, :]
        # x = (x - x[:, 0:20, :].mean(axis=1)) / x[:, 0:20, :].mean(axis=1)
        x *= 100
    return x

#
#  MAIN PLOTTING SCRIPT
#


#  FIGURE 1
if 1 in to_plot:
    X = [np.asarray(smokers[k])[:1000, start_time:] for k in keys]
    X = preprocess_array(X)
    X /= n_smokers

    times = np.arange(X.shape[1])

    fig = plt.figure()
    for i, (key, color) in enumerate(zip(keys, colors_context)):
        sns.tsplot(X[::decim, :, i:i + 1], color=color,
                   condition=mapping[key])

    plt.xlabel('Time [AU]')
    plt.xticks(n_xticks)
    plt.xlim(0, max(n_xticks))
    plt.ylabel('Proportion of smokers')
    if set_titles:
        plt.title('Evolution of Smoking prevalence')
    sns.despine(offset=10, trim=True)
    plt.gca().set_xticklabels(np.linspace(0, 1, len(n_xticks)))
    plt.subplots_adjust(bottom=0.15)
    plt.show()
    fig.savefig(op.join('figures', 'n_smokers-%i.png' % ci), dpi=300)

####################
#  FIGURE 2
if 2 in to_plot:

    centrality = data['centrality']
    centrality_smoker = np.array([
        np.asarray(centrality['smoker'][k])[:1000, start_time:] for k in keys])
    centrality_nosmoker = np.array([
        np.asarray(centrality['non_smoker'][k])[:1000, start_time:]
        for k in keys])
    centrality_smoker = rescale(preprocess_array(centrality_smoker), axis=1)
    centrality_nosmoker = rescale(preprocess_array(centrality_nosmoker),
                                  axis=1)
    centrality_diff = centrality_smoker - centrality_nosmoker

    fig, axes = plt.subplots(1, 2, figsize=(12, 8), sharex=True, sharey=True)
    ax1, ax2 = axes.ravel()
    ax1.set_xlim(0, max(n_xticks))
    ax2.set_xlim(0, max(n_xticks))
    ax1.set_ylim(-40, 30)
    ax2.set_ylim(-40, 30)
    ax1.text(20, 27, 'a', weight='bold', fontsize=24)
    ax2.text(20, 27, 'b', weight='bold', fontsize=24)

    X = centrality_smoker
    for i, (key, color) in enumerate(zip(keys, colors_context)):
        sns.tsplot(X[::decim, :, i:i + 1],
                   color=color, condition=mapping[key], ax=ax1, ci=[95, ci],
                   legend=False)
    if set_titles:
        fig.suptitle('Evolution of eigenvector centrality', fontsize=32)

    plt.xticks(n_xticks)
    # if set_titles:
    ax1.set_title('smokers')
    ax1.set_xlabel('Time [AU]')
    ax1.set_ylabel('Relative change in centrality [percent]')
    ax1.set_xlim(0, max(n_xticks))

    X = centrality_nosmoker
    for i, (key, color) in enumerate(zip(keys, colors_context)):
        sns.tsplot(X[::decim, :, i:i + 1],
                   color=color, condition=mapping[key], ax=ax2, ci=[95, ci],
                   legend=True)

    plt.xticks(n_xticks)
    ax2.set_xticklabels(np.linspace(0, 1, len(n_xticks)))
    # if set_titles:
    ax2.set_title('non-smokers')
    ax2.set_xlabel('Time [AU]')
    ax2.set_xlim(0, max(n_xticks))

    sns.despine(offset=10, trim=True)
    plt.subplots_adjust(bottom=0.15, top=.85, right=.96, left=0.14, wspace=.27)
    plt.show()
    fig.savefig(op.join('figures', 'centrality_smokers--%i.png' % ci), dpi=300)

####################
#  FIGURE 3

X = np.array([np.asarray(data['conditional_prob'][k])[:1000, start_time:, :]
              for k in keys])
X = np.transpose(X, (-1, 1, 2, 0))
X = moving_average(X, ma_window_size, axis=2)

if 3 in to_plot:
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), sharex=False,
                             sharey=False)
    axes = axes.ravel()
    n_yticks = [0, 200, 400, 700]
    for i, (key, color) in enumerate(zip(keys, colors_context)):
        ax = axes[i]
        ax.set_ylim(-100, 700)
        ax.text(50, 710, 'abcd'[i], weight='bold', fontsize=24)
        for degree, (x, subcolor) in enumerate(
                zip(X,
                    sns.light_palette(color, len(X))[::-1]), 1):
            sns.tsplot(rescale(x[::decim, :, i:i + 1] + 1, axis=1),
                       color=subcolor, condition='degree %i' % degree,
                       ax=ax, ci=ci)
        # plt.xticks(times, times[::1000])
        if i == 0:
            ax.set_ylabel('Relative change in CP [percent]', labelpad=10)
        ax.set_title(mapping[key])
        ax.set_xlim(0, max(n_xticks))
        ax.set_xticks(n_xticks)
        if i in (1, 3):
            ax.set_ylabel('')
            ax.set_yticks(n_yticks, ['', '', '', ''])
        if i in (0, 1):
            ax.set_xlabel('')
            ax.set_yticks(n_yticks, ['', '', '', ''])
        if i in (2, 3):
            ax.set_xlabel('Time [AU]')
        ax.set_xticklabels(np.linspace(0, 1, len(n_xticks)))

    if set_titles:
        fig.suptitle('Conditional probability of smoking', fontsize=32)
    sns.despine(offset=10, trim=True)
    plt.subplots_adjust(left=0.1, bottom=0.15, right=0.97, top=0.86,
                        wspace=.41, hspace=.41)
    fig.show()
    fig.savefig(op.join('figures', 'conditional_prob_all-%i.png' % ci),
                dpi=300)
# fig.savefig(op.join('figures', 'conditional_prob_all.png'))
