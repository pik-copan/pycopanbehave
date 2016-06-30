#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
pycopanbehave -- An adaptive network mode of behaviour selection in Python

Copyright (C) 2011--2016 Potsdam Institute for Climate Impact Research
Authors: Jonathan F. Donges <donges@pik-potsdam.de>,
         Carl-Friedrich Schleussner <schleussner@pik-potsdam.de>,
         Denis Engemann <denis.engemann@gmail.com>
URL:     <http://www.pik-potsdam.de/copan/software>

DEFINES DYNAMIC MODEL CORE CLASS CopanBehaveModel
"""

#
#  Imports
#

import numpy as np

#  Import Network class from pyunicorn
from pyunicorn import Network


class CopanBehaveModel(object):
    """
    Implements COPAN:BEHAVE: An adaptive network mode of behaviour selection.

    The CopanBehaveModel class encapsulates the dynamic core of the
    COPAN:BEHAVE model. A detailed mathematical description of the model can
    be found in the following publication:

      - C.-F. Schleussner(#), J.F. Donges(#), D.A. Engemann(#), and
          A. Levermann,
        Co-evolutionary behaviour selection in adaptive social networks
        predicts clustered marginalization of minorities,
        Preprint: arxiv.org:1512.05013 [physics.soc-ph] (2015).
        (#) The first three authors share the lead authorship.

    Parameters
    ----------
    background_proximity_matrix : np.array[N,N] of dtype float
        The background proximity matrix describing social background structure
    agent_characteristics : np.array[N] of dtype float
        Vector of initial individual endogeneous agent characteristics (e.g.
        smoking behaviour)
    agent_properties : np.array[N] of dtype float
        Vector of initial exogeneous individual agent properties (e.g. smoking
        disposition)
    initial_contact_network : pyunicorn.Network object
        Initial contact network
    interaction_probability_function : function
        Function for computing interaction probability matrix from social
        distance matrix
    L : bunch dictionary
        Dictionary containing model parameters
    coupling_instance : string
        Flag determines which model type (full or partial) is used,
        possible values:
        'full' : fully coupled model
        'infl' : local social influence, static contact network
        'mean_field' : mean-field social influence, dynamic contact network
        'dyn' : no social influence, dynamic contact network
    degree_preference : np.array[N] of dtype int
        Vector of individual degree preferences
    """

    #
    #  Class methods
    #

    def __init__(self, background_proximity_matrix, agent_characteristics,
                 agent_properties, initial_contact_network,
                 interaction_probability_function, L, coupling_instance='full',
                 degree_preference=None):
        """
        Constructor of CopanBehaveModel class.

        See class docstring for parameter description.
        """
        #
        #  Initialize instance variables
        #

        #  Background proximity matrix
        self._background_proximity_matrix = background_proximity_matrix.copy()

        #  Vector of initial endogenous individual characteristics
        self._agent_characteristics = agent_characteristics.copy()

        #  Vector of initial exogenous individual properties
        self._agent_properties = agent_properties.copy()

        #  Flag determines which model type (full or partial) is used
        self._coupling_instance = coupling_instance

        #  Initial contact network
        adjacency = initial_contact_network.adjacency
        self._contact_network = Network(adjacency=adjacency, directed=False,
                                        silence_level=3)

        #  Number of interactions (edges in interaction network)
        #  in a model time step
        self._no_interactions = 0.

        #  Number of new edges in contact network appearing
        #  in a model time step
        self._new_edges = 0

        #  Bunch dictionary containing model parameters
        self._L = L

        #  Weight parameters for computation of proximity matrix
        self._char_weight = L.char_weight

        #  Vector of individual degree preferences
        if degree_preference is not None:
            self._degree_preference = degree_preference.copy()
        else:
            #  If no degree preference passed to __init__, use degree vector
            #  of initial contact network
            self._degree_preference = initial_contact_network.degree()

        #  Set interaction_probability function
        self._interaction_probability_function = \
            interaction_probability_function

        #  Set number of agents N
        self._N = self._contact_network.N

    #
    #  Service methods
    #

    def get_background_proximity_matrix(self):
        return self._background_proximity_matrix

    def get_agent_characteristics(self):
        return self._agent_characteristics

    def set_agent_characteristics(self, new_chd):
        self._agent_characteristics = new_chd

    def get_agent_properties(self):
        return self._agent_properties

    def set_agent_properties(self, new_chd):
        self._agent_properties = new_chd

    def get_contact_network(self):
        return self._contact_network

    def get_degree_preference(self):
        return self._degree_preference

    #
    #  Model methods
    #

    def get_proximity_matrix(self, agent_characteristics, char_weight,
                             background_proximity_matrix):
        """
        Returns a proximity matrix based on distance between individual
        characteristics and a background proximity structure.

        Parameters
        ----------
        agent_characteristics : np.array[N] of dtype float
            Vector of initial individual endogeneous agent characteristics
            (e.g. smoking behaviour)
        char_weight : np.array[2]
            Weight parameters for computation of proximity matrix
        background_proximity_matrix : np.array[N,N] of dtype float
            The background proximity matrix describing social background structure

        Returns
        -------
        Proximity matrix : np.array[N,N] of dtype float
        """
        distances = np.zeros((self._N, self._N))

        for i in range(self._N):
            #  Calculate the distances based on the first characteristic
            distances[i, :] = (char_weight[0] * np.abs(
                                   np.abs(agent_characteristics -
                                          agent_characteristics[i]) - 1))

        #  Compute proximity matrix as weighted sum of distance in
        #  characteristics and background proximity structure
        return distances + char_weight[1] * background_proximity_matrix

    def iterate(self, n_steps):
        """
        Iterate the model for n_steps further time steps.

        Parameters
        ----------
        n_steps : int
            Number of model time steps to be integrated
        """
        # Get current contact network
        contact_network = self.get_contact_network()

        # Get degree preference
        degree_preference = self.get_degree_preference()

        # Get interaction_probability function
        interaction_probability_function = \
            self._interaction_probability_function

        # Iterate
        for i in range(n_steps):
            #  Calculate distance metric from contact network
            distance_metric = self.get_distance_metric_matrix(contact_network)

            #  Calculate interaction_probability matrix from distance metric
            #  matrix
            interaction_probability = self.get_interaction_probability_matrix(
                                           distance_metric,
                                           interaction_probability_function)

            #  Obtain interaction network by randomly drawing
            #  using interaction probabilities
            interaction_network = self.get_interaction_network(
                interaction_probability)

            # Implementation of the social influence scheme
            if self._coupling_instance in ['full', 'infl', 'mean_field']:
                agent_characteristics = self.update_social_influence(
                    self._agent_characteristics, self._agent_properties,
                    interaction_network)
            else:
                agent_characteristics = self._agent_characteristics

            # Proximity matrix derived via linear combination of individual
            # characteristics and background proximity structure
            proximity_matrix = self.get_proximity_matrix(
                agent_characteristics,
                self._char_weight,
                self._background_proximity_matrix)

            #  Update contact network
            if self._coupling_instance in ['full', 'dyn', 'mean_field']:
                contact_network = self.update_contact_network(
                    contact_network, proximity_matrix, interaction_network,
                    degree_preference)

        #  Update instance variables to iterated values
        self._agent_characteristics = agent_characteristics
        self._contact_network = contact_network
        self._no_interactions = np.sum(interaction_network)

    def set_eq_network(self):
        """
        Set the equilibrium contact network for present proximity matrix by
        letting everyone interact with everyone else.
        """
        #  Get current contact network
        contact_network = self.get_contact_network()

        #  Get degree preference
        degree_preference = self.get_degree_preference()

        #  Get individual characteristics distribution
        agent_characteristics = self._agent_characteristics

        #  Everyone interacts with everyone else
        interaction_probability_matrix = np.ones((self._N, self._N))

        #  Get proximity matrix from linear combination of distance in
        #  individual characteristics and background proximity structure
        proximity = self.get_proximity_matrix(
            agent_characteristics, self._char_weight,
            self._background_proximity_matrix)

        #  Update contact network
        contact_network = self.update_contact_network(
            contact_network, proximity, interaction_probability_matrix,
            degree_preference)

        #  Update instance variables to iterated values
        self._contact_network = contact_network
        self._no_interactions = np.sum(interaction_probability_matrix)

    def get_distance_metric_matrix(self, contact_network):
        """
        Returns the geodesic graph distance (shortest path length) matrix
        as a measure of social distance on the contact network.

        Parameters
        ----------
        contact_network : pyunicorn.Network object
            Contact network

        Returns
        -------
        Distance metric on contact network : np.array[N,N] of dtype float
        """
        #  Return shortest path lengths matrix
        return contact_network.path_lengths()

    def get_interaction_probability_matrix(self, distance_metric_matrix,
                                           interaction_probability_function):
        """
        Returns the interaction probability matrix given a distance metric on
        the contact network using the given interaction_probability function.

        Parameters
        ----------
        distance_metric_matrix : np.array[N,N] of dtype float
            Distance metric on contact network
        interaction_probability_function : function
            Function for computing interaction probability matrix from
            social distance matrix

        Returns
        -------
        Interaction probability matrix : np.array[N,N] of dtype float
                                            in interval [0,1]
        """
        return interaction_probability_function(distance_metric_matrix,
                                                self._L)

    def get_interaction_network(self, interaction_probability_matrix):
        """
        Returns the interaction network's adjacency matrix stochstically given
        the interaction probability matrix.

        Parameters
        ----------
        interaction_probability_matrix : np.array[N,N] of dtype float
            Interaction probability matrix (in interval [0,1])

        Returns
        -------
        Interaction network's adjacency matrix : np.array[N,N] of dtype int
        """
        #  Draw uniformly distributed random numbers from the interval [0,1]
        random_numbers = np.random.rand(self._N, self._N)
        #  Symmetrize
        random_numbers = (random_numbers + random_numbers.transpose()) / 2.

        #  Return adjacency matrix of interaction_probability network
        return (random_numbers <=
                interaction_probability_matrix).astype("int8")

    def update_social_influence(self, agent_characteristics, agent_properties,
                                interaction_network):
        """
        Return updated agent characteristics following an Ising-type scheme.

        Parameters
        ----------
        agent_characteristics : np.array[N] of dtype float
            Vector of initial individual endogeneous agent characteristics
            (e.g.smoking behaviour)
        agent_properties : np.array[N] of dtype float
            Vector of initial exogeneous individual agent properties
            (e.g. smoking disposition)
        interaction_network : np.array[N,N] of dtype int
            Interaction network's adjacency matrix

        Returns
        -------
        Updated individual agent characteristics : np.array[N] of float
        """
        #  Initialize
        agent_characteristics_update = agent_characteristics

        #  Compute share of smokers in the system
        sm_share = float(agent_characteristics.sum()) / self._N

        for i in range(len(agent_characteristics)):
            #  Get number of interactions of node i
            nai = np.sum(interaction_network[i, :])

            #  No update takes place when there are no interactions
            if nai == 0:
                agent_characteristics_update[i] = agent_characteristics[i]
            #  Update of individual characteristics
            else:
                random_number = np.random.rand(1)
                #  Update of characteristics is based on social influence by
                #  the actual peers in the interaction network
                if self._coupling_instance in ['full', 'infl']:
                    #  Agent i is non-smoker
                    if agent_characteristics[i] == 0:
                        prob_flip = (self._L.C * agent_properties[i] *
                                     np.sum(interaction_network[i, :] *
                                            agent_characteristics) / nai)
                        agent_characteristics_update[i] = (
                            random_number <= prob_flip).astype("int8")
                    #  Agent i is smoker
                    else:
                        prob_flip = (self._L.C * (1 - agent_properties[i]) *
                                     (1 - np.sum(interaction_network[i, :] *
                                      agent_characteristics) / nai))
                        agent_characteristics_update[i] = \
                            1 - (random_number <= prob_flip).astype("int8")

                #  Mean field social influence for the dyn only case
                elif self._coupling_instance in ['mean_field']:
                    #  Agent i is non-smoker
                    if agent_characteristics[i] == 0:
                        prob_flip = self._L.C*agent_properties[i] * sm_share
                        agent_characteristics_update[i] = (
                            random_number <= prob_flip).astype("int8")
                    #  Agent i is smoker
                    else:
                        prob_flip = (self._L.C *
                                     (1 - agent_properties[i]) *
                                     (1 - sm_share))
                        agent_characteristics_update[i] = \
                            1 - (random_number <= prob_flip).astype("int8")

        return agent_characteristics_update

    def update_contact_network(self, contact_network, proximity_matrix,
                               interaction_network, degree_preference):
        """
        Returns the updated contact network.

        Parameters
        ----------
        contact_network : pyunicorn.Network object
            Contact network
        proximity_matrix : np.array[N,N] of dtype float
            Proximity matrix
        interaction_network : np.array[N,N] of dtype int
            Interaction network's adjacency matrix
        degree_preference : np.array[N] of dtype int
            Vector of individual degree preferences

        Returns
        -------
        Updated contact network : pyunicorn.Network object
        """
        #  Get old contact adjacency matrix
        old_contact_adjacency = contact_network.adjacency

        #  Initialize new contact network adjacency matrix
        contact_adjacency = np.zeros(old_contact_adjacency.shape, dtype="int8")

        # Initialize list of sorted indices
        potential_contact_indices = []

        # Loop over all agents
        for i in range(self._N):
            # Get node indices of contacts
            contact_indices = np.where(old_contact_adjacency[i, :] == 1)[0]

            # Agent i has contacts
            # if k[i] != 0:  Get node indices of interaction network
            interaction_probability_indices = np.where(
                interaction_network[i, :] == 1)[0]

            # Combine both lists of indices and discard repeated entries
            indices = list(np.unique(np.append(contact_indices,
                                     interaction_probability_indices)))

            # Get proximity values of contacts and interaction_network
            similarities = proximity_matrix[i, indices]

            indices = np.array(indices)

            # Get degree_preference[i] indices with largest similarities
            sorted_indices = indices[
                similarities.argsort()][-degree_preference[i]:]

            #  Get k indices with largest similarities
            potential_contact_indices.append(sorted_indices)

            # #  Agent i has no contacts
            # else:   potential_contact_indices.append([])

        #  Keep only bidirectional potential new contacts
        for i in range(self._N):
            filtered_indices = []

            for j in potential_contact_indices[i]:
                if i in potential_contact_indices[j]:
                    filtered_indices.append(j)
            contact_adjacency[i, filtered_indices] = 1

        # return the no of changed changed contacts per turn (
        # nodes can have edges with themselves, the diagonal is therefore
        # removed)
        self._new_edges = (
            (np.sum(np.abs(old_contact_adjacency - contact_adjacency)) -
             np.sum(np.diag(np.abs(old_contact_adjacency -
                                   contact_adjacency)))) / 2)

        #  Update contact Network object
        contact_network.adjacency = contact_adjacency

        #  Return new contact network
        return contact_network
