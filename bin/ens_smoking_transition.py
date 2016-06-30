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

MAIN FUNCTIONS FOR SMOKING BEHAVIOUR TRANSITION SIMULATIONS
"""

#
#  Imports
#

import sys
import pickle

import numpy as np
import scipy.stats
import scipy.integrate as integ

import igraph

# Add repository path of the model core
sys.path.append('../')

# Import pyunicorn dependencies
from pyunicorn import Network

#  Import CopanBehaveModel class from pycopanbehave
from pycopanbehave import CopanBehaveModel


#
#  Define Bunch class
#

class Bunch(dict):
    """Bunch object"""
    def __init__(self,**kw):
        dict.__init__(self,kw)
        self.__dict__.update(kw)


#range = getattr(__builtins__, 'xrange', __builtins__.range)


###############################################################################
#
#  Functions
#
###############################################################################


def distribution_function(yb, x):
    """"
    Modulate the parabolic distribution y= a(b-x)**2 +c fct. with the
    following boundary conditions:
      - 1. x=0: y == 3 : follows from the initial distribution function
      - 2. Int[0,1] == 1: normalization criteria
      - 3. x=1: y == yb, where yb is the tunable parameter, which reflects
          the probability for an agent to have smoking disposition 1
        yb0==3
    """
    # print('yb', yb)
    # yb = 3
    b = (1. + yb / 3) / (1 + yb)
    a = 2. / (b - 1. / 3)
    c = 3. - a * b ** 2

    return a * (b - x) ** 2 + c


def rejection_sampling(N, distribution_function_ini, yb):
    """
    Creates random sampling of an arbitrary distribution using the rejection
    sampling method for a continuous characteristic (e.g. smoking disposition)
    and a second binary characteristic based on this.
    """
    result = np.zeros((2, N))

    i = 0
    while i < N:
        random_numbers = np.random.rand(2)
        if random_numbers[0] < distribution_function_ini(yb,
                                                         random_numbers[1]):
            result[0, i] = random_numbers[1]
            # if the smoking disposition is greater than a certain value,
            # the binary characteristic smoker/non-smoker is assigned
            result[1, i] = (random_numbers[1] > .5).astype('int8')
            i += 1

    return result[0, :], result[1, :]


def interaction_probability_function(distance_metric_matrix, L):
    """
    Returns the interaction probability matrix given a distance metric on the
    contact network.

    The interaction probability is computed applying a exp(-d/9)-relationship
    according to the three degrees of influence findings by ... .

    Since the number of nodes with a certain path length is gaussian
    distributed with a peak around the mean average path length, this has to
    be accounted for. Thereby, we norm the interaction probability according
    to the absolute number of occurences of each entry.
    """
    # Three degrees of influence plus an interaction offset to avoid
    # disconnection
    exp_dec = (L.p_ai-L.interaction_offset) * \
        np.exp(-(distance_metric_matrix-1) / 2.) + L.interaction_offset
    distmax = distance_metric_matrix[np.isfinite(distance_metric_matrix)].max()
    histo_range = np.arange(1, distmax)
    distribution = np.histogram(distance_metric_matrix.flatten(), histo_range)

    for i in distribution[1][:-1]:
        exp_dec[distance_metric_matrix == i] *= (float(distribution[0][0]) /
                                                 distribution[0][i-1])

    # exp_dec=np.exp(-(distance_metric_matrix-1)/2.)
    return exp_dec


def transient_disposition_distribution(N, disp_distr_t, yb):
    """
    This function governs the transition between two distributions.

    The functions' similarity is derived via a Kolmogorov-Smirnov test
    as a necessary criterium
    (https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test).

    The target is set to 0.1, which is equivalent to the value one gets
    for random sampling from the given distribution.

    The noise added is lognormal distributed to allow for
    long-range jumps and scaled with the deviation of the actual distribution
    to the ideal curve. The greater the positive deviation, the greater the
    noise.
    """

    #  Define helper function
    def integrate_cdf(input_array):
        cdf = np.zeros(input_array.shape[0])
        for i in range(input_array.shape[0]):
            cdf[i] = integ.quad(lambda x: distribution_function(yb, x),
                                0, input_array[i])[0]
        return cdf

    kolm_smir_step = 2

    # Kolm Smir target
    target = 0.1
    # noise scaling coeff
    tr_noise_coeff = 0.1

    # set counters
    k = 1
    n_inc = 0

    # derive the ideal distribution for given yb
    x = np.linspace(0, 1, 101)
    target_dist = distribution_function(yb, x)
    target_dist = target_dist * float(N) / (sum(target_dist))
    # get deviation
    hist_values_sample, hist_bins = np.histogram(
        disp_distr_t, bins=100, range=(0, 1))
    hist_values_sample = np.append(hist_values_sample, hist_values_sample[-1])
    hist_diff = hist_values_sample - target_dist
    # get the noise level for all N agents
    disp_distr_round = np.asarray(100 * disp_distr_t, dtype='int')
    distr_dev = hist_diff[disp_distr_round]
    distr_dev = np.asarray(distr_dev, dtype='float')
    # scale and set positive
    distr_dev = distr_dev / 50
    distr_dev += np.abs(distr_dev.min())
    logn_std = 1.
    while kolm_smir_step >= target:
        # generate random lognormal dist and sign
        random_noise_increase = np.random.lognormal(0, logn_std, N)
        random_sign = np.random.randint(-1, 1, N)
        random_sign[random_sign == 0] = 1
        # add the noise
        disp_distr_tp1 = (disp_distr_t + tr_noise_coeff * random_sign *
                          random_noise_increase * distr_dev)
        # check for boundaries
        disp_distr_tp1[disp_distr_tp1 > 1] = disp_distr_t[disp_distr_tp1 > 1]
        disp_distr_tp1[disp_distr_tp1 < 0] = disp_distr_t[disp_distr_tp1 < 0]
        # Kolmogorov-Smirnov test
        kolm_smir_step = scipy.stats.kstest(disp_distr_tp1, integrate_cdf)[0]
        # print('kolm_smir_step',kolm_smir_step)
        k += 1
        # in case of non-convergence, increase the standard deviation for the
        # lognormal dist to allow for
        if k > 100:
            logn_std = logn_std * 1.2
            print('increased noise, logn_std', logn_std)
            k = 1
            n_inc += 1
            if n_inc > 10:
                print('DISTRIBUTION TRANSITION FAILED',
                      'kolm_smir_step', kolm_smir_step, 'yb', yb)
                return disp_distr_t
    return disp_distr_tp1


def generate_initial_distance_sm(L):
    """
    The initial proximity structure is generated. It is a small world NW
    following the classical Watts & Strogatz approach with an initial ring
    that is successively rewired with the rewiring probability L.p_rew.
    According to Watts & Strogatz 0.01 < p_rew < 0.1.
    """
    proximity_small_world = np.zeros((L.N, L.N))
    proxim_nw = igraph.GraphBase.Lattice([L.N], nei=L.links_per_node,
                                         directed=False, mutual=True,
                                         circular=True)
    proxim_nw.rewire(int(L.p_rew * L.N * L.links_per_node))
    small_world_distance_matrix = np.asarray(proxim_nw.shortest_paths())

    # Derive the proximity matrix on the basis of the small-world network,
    # random noise added to generate a continuum
    random_prox_noise = np.random.random((L.N, L.N))
    random_prox_noise = (random_prox_noise+random_prox_noise.transpose()) * .5
    proximity_small_world = (1 - .1 * (small_world_distance_matrix-1) - 0.1 *
                             random_prox_noise)

    # introduce a proximity offset,
    # if the distance falls below 5, there should be no considerable difference
    # anymore
    proximity_small_world[np.where(proximity_small_world <= 0.2)] = 0.2

    # ensure that the diagonals are 0
    for k in range(proximity_small_world.shape[0]):
        proximity_small_world[k, k] = 0

    return proximity_small_world


def calc_cond_prob(smokers, nw_full, deg_sep_max, N):
    """
    Add docstring!
    """
    rcp = np.zeros(5)
    for i in range(deg_sep_max):
        deg_sep = i + 1
        smoking_dep = []
        for node in smokers:
            distance_matrix = nw_full.path_lengths()
            contact_one = np.where(distance_matrix[node, :] == deg_sep)
            if contact_one[0].size > 0:
                smoking_dep.append(
                    np.sum(nw_full.node_attribute('smoker')[contact_one]) /
                    float(contact_one[0].size) / (float(len(smokers)) / N) - 1)
        rcp[i] = np.mean(smoking_dep)
    return rcp


def derive_nw_chars(outdic, model_trans, L, i=0):
    """
    Add docstring!
    """
    if i == 0:
        l = L.n_iterations
        # Initialize list of full output variables implemented
        nw_chars_keylist = [
            'no_interactions',
            'cogn_diss_coeff',
            'smokers',
            'no_of_smokers',
            'mean_degree',
            'clustering',
            'apl',
            'contact_change',
            'no_disconnected_nodes',
            'kolm_smir',
            'conditional_prob',
            'disp_distr',
            'max_sm_cl_size',
            'ave_sm_cl_size']
        sm_chars_list = [
            'centrality',
            'betweenness',
            'closeness',
            'degree',
            'ind_clustering']
        other_vars_list = ['conditional_prob', 'disp_distr']

        # Check, if full output should written or just selection
        if L.write_full_output_to_file:
            nw_chars_keylist_out = nw_chars_keylist
            sm_chars_list_out = sm_chars_list
            other_vars_list_out = other_vars_list
        else:
            nw_chars_keylist_out = [keys for keys in nw_chars_keylist if
                                    keys in L.variables_out]
            sm_chars_list_out = [keys for keys in sm_chars_list if keys in
                                 L.variables_out]
            other_vars_list_out = [keys for keys in other_vars_list if keys in
                                   L.variables_out]

        ######################
        # initialize output dictionary
        ######################
        for k in nw_chars_keylist_out:
            outdic[k] = np.zeros(l)

        for k in sm_chars_list_out:
            outdic[k] = {'smoker': np.zeros(l), 'non_smoker': np.zeros(l)}

        if 'conditional_prob' in other_vars_list_out:
            outdic['conditional_prob'] = np.zeros((l, 5))
        if 'disp_distr' in other_vars_list_out:
            outdic['disp_distr'] = np.zeros((l, L.N))

    ###########################################################################
    # write output model chars
    smokers = np.where(model_trans.get_agent_characteristics() == 1)[0]

    non_smokers = np.where(model_trans.get_agent_characteristics() == 0)[0]
    if 'no_interactions' in outdic:
        outdic['no_interactions'][i] = model_trans._no_interactions
    if 'contact_change' in outdic:
        outdic['contact_change'][i] = model_trans._new_edges
    if 'no_of_smokers' in outdic:
        outdic['no_of_smokers'][i] = len(smokers)
    if 'cogn_diss_coeff' in outdic:
        outdic['cogn_diss_coeff'][i] = np.mean(
            np.abs(model_trans.get_agent_characteristics() -
                   model_trans.get_agent_properties()))
    if 'disp_distr' in outdic:
        outdic['disp_distr'][i, :] = model_trans.get_agent_properties()

    # write output network chars
    if 'mean_degree' in outdic:
        outdic['mean_degree'][i] = (model_trans.get_contact_network()
                                               .degree()
                                               .mean())
    if 'clustering' in outdic:
        outdic['clustering'][i] = (model_trans.get_contact_network()
                                              .global_clustering())

    if 'apl' in outdic:
        outdic['apl'][i] = (model_trans.get_contact_network()
                                       .average_path_length())

    if 'ave_sm_cl_size' in outdic:
        #######
        # derive max smoker cluster size
        ts_adj = model_trans.get_contact_network().adjacency
        sm_mat = np.zeros(1000)
        sm_mat[smokers] = 1
        smok = np.asarray(sm_mat, dtype='bool')

        try:
            nwsm = Network(adjacency=ts_adj[smok, :][:, smok])
            get_cl = nwsm.graph.clusters(mode='WEAK')
            sm_cl_len = [len(get_cl[rft]) for rft in range(len(get_cl))]
            sm_cl_len = np.asarray(sm_cl_len)
            outdic['max_sm_cl_size'][i] = np.max(sm_cl_len)
            outdic['ave_sm_cl_size'][i] = np.mean(sm_cl_len[sm_cl_len > 1])
        except:
            print('derive smoker_cluster_size failed')
    if 'centrality' in outdic:
        outdic['conditional_prob'][i, :] = calc_cond_prob(
            smokers, model_trans.get_contact_network(), L.cond_prob_degree,
            L.N)

    if 'centrality' in outdic:
        outdic['centrality']['smoker'][i] = np.mean(
            np.asarray(model_trans.get_contact_network()
                                  .graph
                                  .evcent(scale=False))[smokers])
        outdic['centrality']['non_smoker'][i] = np.mean(
            np.asarray(model_trans.get_contact_network()
                                  .graph.evcent(scale=False))[non_smokers])

    if 'betweenness' in outdic:
        outdic['betweenness']['smoker'][i] = np.mean(
            model_trans.get_contact_network()
                       .betweenness()[smokers])
        outdic['betweenness']['non_smoker'][i] = np.mean(
            model_trans.get_contact_network()
                       .betweenness()[non_smokers])
    if 'closeness' in outdic:
        outdic['closeness']['smoker'][i] = np.mean(
            model_trans.get_contact_network()
                       .closeness()[smokers])
        outdic['closeness']['non_smoker'][i] = np.mean(
            model_trans.get_contact_network()
                       .closeness()[non_smokers])
    if 'degree' in outdic:
        degree_matrix = model_trans.get_contact_network().degree()
        outdic['degree']['smoker'][i] = np.mean(degree_matrix[smokers])
        outdic['degree']['non_smoker'][i] = np.mean(degree_matrix[non_smokers])
    if 'ind_clustering' in outdic:
        outdic['ind_clustering']['smoker'][i] = np.mean(
            model_trans.get_contact_network()
                       .local_clustering()[smokers])
        outdic['ind_clustering']['non_smoker'][i] = np.mean(
            model_trans.get_contact_network()
                       .local_clustering()[non_smokers])
    if 'no_disconnected_nodes' in outdic:
        disc_node = np.where(
            np.isinf(model_trans.get_contact_network()
                                .path_lengths()[1, :]) == True)[0]
        outdic['no_disconnected_nodes'][i] = disc_node.shape[0]

    return outdic

###############################################################################
# SINGLE RUN
###############################################################################


def generate_eq(L):
    """
    Add docstring!
    """
    print('############## INITIALIZE MODEL')
    no_components = 2
    while no_components != 1:

        agent_properties, agent_characteristics = rejection_sampling(
            L.N, distribution_function, L.yb_initial)
        # Degree preference distribution
        degree_preference = np.random.normal(
            L.mean_degree_pref, L.std_degree_pref, L.N).astype("int8")

        # Get the underlying network structure
        background_proximity = generate_initial_distance_sm(L)

        ###################
        # Generate equilibrium contact networks
        ###################

        #  Create random contact network (Erdoess-Renyi graph)
        contact_network = Network.ErdosRenyi(
            n_nodes=L.N, link_probability=L.mean_degree_pref / float(L.N - 1))

        #  Set silence level
        contact_network.silence_level = 3

        # Initialize fully coupled model
        coupling_instance = 'full'

        model_initial = CopanBehaveModel(
            background_proximity, agent_characteristics, agent_properties,
            contact_network, interaction_probability_function, L,
            coupling_instance, degree_preference)

        model_initial.set_eq_network()
        (model_initial.get_contact_network()
                      .set_node_attribute(
                          'smoker',
                          model_initial.get_agent_characteristics()))

        print('############## RUN EQ GENERATION')
        print('EQ generation', sum(model_initial.get_agent_characteristics()))
        for i in range(L.n_initial_eq_it):
            model_initial.iterate(1)
            (model_initial.get_contact_network()
                          .set_node_attribute(
                              'smoker',
                              model_initial.get_agent_characteristics()))
        print('--------------- EQ Network generated --------------- ')

        clusters = (model_initial.get_contact_network()
                                 .graph
                                 .clusters(mode='WEAK'))
        no_components = len(clusters[0]) + len(clusters) - L.N
        print('################# Number of components initial network '
              '#################')
        print(no_components, len(clusters))
    return model_initial


def transition(coupling_instance, model_initial, L, out, char_dist):
    """
    Add docstring!
    """
    model_trans = CopanBehaveModel(
        model_initial.get_background_proximity_matrix(),
        model_initial.get_agent_characteristics(),
        model_initial.get_agent_properties(),
        model_initial.get_contact_network(),
        interaction_probability_function,
        L, coupling_instance, model_initial._degree_preference)

    ######################
    # Derive the final distribution
    ######################

    print('transient part')
    (model_trans.get_contact_network()
                .set_node_attribute(
                    'smoker',
                    model_trans.get_agent_characteristics()))
    output_dict = {}
    nw_snapshots_dic = {}
    for i in range(L.n_iterations):

        ###############################################################
        #    TRANSIENT CHANGE OF THE SMOKING DISPOSITION
        ###############################################################

        model_trans.set_agent_properties(char_dist[i, :])
        if coupling_instance == 'dyn':
            model_trans.set_agent_characteristics(
                np.round(model_trans.get_agent_properties(), 0))

        ##############################################################
        model_trans.iterate(1)
        (model_trans.get_contact_network()
                    .set_node_attribute(
                        'smoker',
                        model_trans.get_agent_characteristics().copy()))
        output_dict = derive_nw_chars(output_dict, model_trans, L, i)

        #######################################################################
        # SAVE NETWORK INSTANCES
        #######################################################################

        if ((coupling_instance == 'full') and
            (i == L.n_iterations - 1) or
            (i in np.arange(0, L.n_iterations,
                            L.nw_save_steps, dtype='int')) and
                L.write_full_output_to_file):
            nw_snapshots_dic[i] = (model_trans.get_contact_network()
                                              .adjacency)

    if (coupling_instance == 'full') and (L.write_full_output_to_file):
        with open(L.output_path + '/nw_snaps_eq_yb_%s_%s.pkl' % (
                  L.yb_final, out), 'wb') as fid:
            pickle.dump(nw_snapshots_dic, fid)

    return output_dict


def do_one(out, L):
    """
    Add docstring!
    """
    print('out', out)

    ###########################################################################
    #
    #  MAIN SCRIPT
    #
    ###########################################################################

    # INITIALIZE NETWORK INSTANCE
    model_initial = generate_eq(L)
    out_dic_2x2 = {}

    ###############################################################
    #    TRANSIENT CHANGE OF THE SMOKING DISPOSITION
    ###############################################################
    char_dist = np.empty(((L.n_iterations, L.N)))
    char_dist[0, :] = model_initial.get_agent_properties()
    for i in range(1, L.n_transition):
        yb = (L.yb_initial - (L.yb_initial - L.yb_final) * (float(i + 1) /
              L.n_transition))
        char_dist[i, :] = transient_disposition_distribution(
            L.N, char_dist[i-1, :], yb)
    if L.write_char_disposition:
        fname = L.output_path + '/chararacter_disposition_%s.pkl' % out
        with open(fname, 'wb') as fid:
            pickle.dump(char_dist, fid)
    print('############## SMOKING DISPOSITION MATRIX SUCCESSFULLY CREATED')

    if L.coupling_instances['full']:
        print('############## RUN FULLY COUPLED')
        out_dic_2x2['full'] = transition(
            'full', model_initial, L, out, char_dist)
    if L.coupling_instances['infl']:
        print('############## RUN SOCIAL INFLUENCE ONLY')
        out_dic_2x2['infl'] = transition(
            'infl', model_initial, L, out, char_dist)
    if L.coupling_instances['dyn']:
        print('############## RUN DYNAMIC ONLY')
        out_dic_2x2['dyn'] = transition(
            'dyn', model_initial, L, out, char_dist)
    if L.coupling_instances['mean_field']:
        print('############## RUN MEAN FIELD FORCING')
        out_dic_2x2['mean_field'] = transition(
            'mean_field', model_initial, L, out, char_dist)
    fname = L.output_path + '/trans_smok_%s.pkl' % out
    with open(fname, 'wb') as fid:
        pickle.dump(out_dic_2x2, fid)
