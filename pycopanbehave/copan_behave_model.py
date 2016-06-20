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
	"""

	#
	#  Class methods
	#

	def __init__(self, background_proximity_matrix, agent_characteristics,agent_properties,
	             initial_contact_network, interaction_probability_function, L,
	             coupling_instance='full', degree_preference=None):
		"""
		Constructor of CopanBehaveModel class.

		Add explanation of class constructor arguments here!
		"""
		#
		#  Initialize instance variables
		#

		# Reference proximity to which the model is relaxed to
		self._background_proximity_matrix = background_proximity_matrix.copy()

		# Distribution of the endogenous individual characteristics 
		self._agent_characteristics = agent_characteristics.copy()
		
		# Distribution of the exogenous individual properties 
		self._agent_properties = agent_properties.copy()

		# Flags if the character feedback is active or not
		self._coupling_instance = coupling_instance

		adjacency = initial_contact_network.adjacency
		self._contact_network = Network(adjacency=adjacency, directed=False,
										silence_level=3)

		self._no_interactions = 0.

		self._new_edges = 0

		self._L = L

		self._char_weight = L.char_weight


		if degree_preference is not None:
		    self._degree_preference = degree_preference.copy()
		else:
		    #  If no degree preference passed to __init__, use degree vector
		    #  of initial contact network
		    self._degree_preference = initial_contact.degree()

		#  Set interaction_probability function
		self._interaction_probability_function = interaction_probability_function

		#  Set number of agents N
		self._N = self._contact_network.N

		self._background_proximity_matrix = background_proximity_matrix

	#
	#  Service methods
	#

	def get_background_proximity_matrix(self):
	    return self._background_proximity_matrix

	def get_agent_characteristics(self):
	   return self._agent_characteristics

	def set_agent_characteristics(self,new_chd,):
	   self._agent_characteristics=new_chd

	def get_agent_properties(self):
	   return self._agent_properties

	def set_agent_properties(self,new_chd):
	   self._agent_properties=new_chd

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
		"""
		distances = np.zeros((self._N, self._N))

		for i in xrange(self._N):
			# Calculate the distances based on the first two chars that apply
			# to all nodes uniformly
			distances[i,:] = char_weight[0] \
				* np.abs(np.abs(agent_characteristics- agent_characteristics[i]) - 1) \
				+ char_weight[1] \
				* np.abs(np.abs(agent_characteristics - agent_characteristics[i]) - 1)

		return distances + char_weight[2] * background_proximity_matrix

	def iterate(self, n_steps):
		"""
		Iterate the model for n_steps further time steps.
		"""
		#  Get current contact network
		contact_network = self.get_contact_network()

		#  Get degree preference
		degree_preference = self.get_degree_preference()

		#  Get interaction_probability function
		interaction_probability_function = self._interaction_probability_function

		#  Iterate
		for i in xrange(n_steps):
			#  Calculate distance metric from contact network
			distance_metric = self.get_distance_metric_matrix(contact_network)

			#  Calculate interaction_probability matrix from distance metric
			#  matrix
			interaction_probability = self.get_interaction_probability_matrix(
										   distance_metric,
										   interaction_probability_function)

			#  Obtain interaction network by randomly drawing
			#  using interaction probabilities
			interaction_network = self.get_interaction_network(interaction_probability)

			# Implementation of the social influence scheme
			if self._coupling_instance in ['full','infl','mean_field']:
				agent_characteristics = self.update_social_influence(
										 self._agent_characteristics,
										 self._agent_properties,
										 interaction_network)
			else:
				agent_characteristics = self._agent_characteristics

			# Proximity matrix derived via linear combination of individual
			# characteristics and background proximity structure
			proximity_matrix = self.get_proximity_matrix(agent_characteristics,
											self._char_weight,
							 				self._background_proximity_matrix)

			#  Update contact network
			if self._coupling_instance in ['full','dyn','mean_field']:
				contact_network = self.update_contact_network(contact_network,
											proximity_matrix,
											interaction_network,
											degree_preference)

		#  Update instance variables to iterated values
		self._agent_characteristics = agent_characteristics
		self._contact_network = contact_network
		self._no_interactions = np.sum(interaction_network)

	def get_eq_network(self):
		"""
		Get the optimal target network for a given proximity matrix by
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
		proximity = self.get_proximity_matrix(agent_characteristics,
											self._char_weight,
											self._background_proximity_matrix)

		#  Update contact network
		contact_network = self.update_contact_network(contact_network,
											proximity,
											interaction_probability_matrix,
											degree_preference)

		#  Update instance variables to iterated values
		self._agent_characteristics = agent_characteristics
		self._contact_network = contact_network
		self._no_interactions = np.sum(interaction_probability_matrix)

	def get_distance_metric_matrix(self, contact_network):
	    """
	    Returns the geodesic graph distance (shortest path length) matrix
	    as a measure of social distance on the contact network.
	    """
	    #  Return shortest path lengths matrix
	    return contact_network.path_lengths()

	def get_interaction_probability_matrix(self, distance_metric_matrix,
										   interaction_probability_function):
	    """
	    Returns the interaction probability matrix given a distance metric on
		the contact network using the given interaction_probability function.
	    """
	    return interaction_probability_function(distance_metric_matrix,
	    										self._L)

	def get_interaction_network(self, interaction_probability_matrix):
	    """
	    Returns the interaction network's adjacency matrix given
	    the interaction probability matrix.
	    """
	    #  Draw uniformly distributed random numbers from the interval [0,1]
	    random_numbers = np.random.rand(self._N, self._N)
	    #  Symmetrize
	    random_numbers = (random_numbers + random_numbers.transpose()) / 2.

	    #  Return adjacency matrix of interaction_probability network
	    return (random_numbers <= interaction_probability_matrix).astype("int8")

	def update_social_influence(self, agent_characteristics,agent_properties,
							 interaction_network):
		"""
		Updated the agent characteristics following an Ising-type update scheme
		"""
		#  Initialize
		agent_characteristics_update = agent_characteristics

		#  Compute share of smokers in the system
		sm_share = float(agent_characteristics.sum()) / self._N


		# print 'agent_characteristics_update ini', sum(agent_characteristics_update)
		
		# ave_prob_flip=[]
		for i in xrange(len(agent_characteristics)):
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
				if self._coupling_instance in ['full','infl']:
					#  Agent i is non-smoker
					if agent_characteristics[i] == 0:
						prob_flip = self._L.C*agent_properties[i]*np.sum(interaction_network[i,:]*agent_characteristics)/nai
						agent_characteristics_update[i] = (random_number <= prob_flip).astype("int8")
					#  Agent i is smoker
					else:
						prob_flip = self._L.C*(1-agent_properties[i])*(1-np.sum(interaction_network[i,:]*agent_characteristics)/nai)
						agent_characteristics_update[i] = 1-(random_number <= prob_flip).astype("int8")
						# ave_prob_flip.append(prob_flip)

				#  Mean field social influence for the dyn only case
				elif self._coupling_instance in ['mean_field']:
					#  Agent i is non-smoker
					if agent_characteristics[i] == 0:
						prob_flip = self._L.C*agent_properties[i]*sm_share
						agent_characteristics_update[i] = (random_number <= prob_flip).astype("int8")
					#  Agent i is smoker
					else:
						prob_flip = self._L.C*(1-agent_properties[i])*(1-sm_share)
						agent_characteristics_update[i] = 1-(random_number <= prob_flip).astype("int8")
						# ave_prob_flip.append(prob_flip)

		# print 'agent_characteristics_update after', sum(agent_characteristics_update)
		# print self._coupling_instance, 'ave_prob_flip ',np.mean(np.asarray(ave_prob_flip)),' share of smokers ', sm_share
		return agent_characteristics_update

	def update_contact_network(self, contact_network, proximity_matrix,
						   interaction_network, degree_preference):
		"""
		Returns the updated contact network.
		"""
		#  Get old contact adjacency matrix
		old_contact_adjacency = contact_network.adjacency

		#  Get old degree k
		k = contact_network.degree()

		#  Initialize new contact network adjacency matrix
		contact_adjacency = np.zeros(old_contact_adjacency.shape, dtype="int8")

		#  Initialize list of sorted indices
		potential_contact_indices = []

		#  Loop over all agents
		for i in xrange(self._N):
			#  Get node indices of contacts
			contact_indices = np.where(old_contact_adjacency[i,:] == 1)[0]

			#  Agent i has contacts
			# if k[i] != 0:
			#  Get node indices of interaction network
			interaction_probability_indices = np.where(interaction_network[i,:] == 1)[0]

			#  Combine both lists of indices and discard repeated entries
			indices = list(np.unique(np.append(contact_indices,
			           		interaction_probability_indices)))

			#  Get proximity values of contacts and interaction_network
			similarities = proximity_matrix[i, indices]

			indices = np.array(indices)

			#  Get degree_preference[i] indices with largest similarities
			sorted_indices = indices[similarities.argsort()][-degree_preference[i]:]

			#  Get k indices with largest similarities
			potential_contact_indices.append(sorted_indices)

			# #  Agent i has no contacts
			# else:
				# potential_contact_indices.append([])

		#  Keep only bidirectional potential new contacts
		for i in xrange(self._N):
			# if k[i] != 0:
			#  Initialize
			filtered_indices = []

			for j in potential_contact_indices[i]:
			    if i in potential_contact_indices[j]:
			        filtered_indices.append(j)

			contact_adjacency[i,filtered_indices] = 1


		# return the no of changed changed contacts per turn (nodes can have edges with themselves, the diagonal is therefore removed)
		self._new_edges = (np.sum(np.abs(old_contact_adjacency-contact_adjacency)) - np.sum(np.diag(np.abs(old_contact_adjacency-contact_adjacency)))) / 2
		#print 'new edges:',np.sum((old_contact_adjacency-contact_adjacency))

		#print "Symmetry:", ((contact_adjacency - contact_adjacency.transpose())**2).sum()

		#  Update contact Network object
		contact_network.adjacency = contact_adjacency

		#  Return new contact network
		return contact_network
