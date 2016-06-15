#!/usr/bin/python
# -*- coding: utf-8 -*-
#
"""
pynamicsociety -- Dynamic social network models in Python
Copyright (C) 2013 Potsdam Institute for Climate Impact Research
Authors: Jonathan F. Donges <donges@pik-potsdam.de>,
         Carl-Friedrich Schleussner <schleussner@pik-potsdam.de>
URL: <http://www.pik-potsdam.de/members/donges/software>
"""

#
#  Imports
#

#  Import numpy for numerics
import numpy as np

#  Import Network class from pyunicorn
from pyunicorn import Network


class DynamicSocietyModel(object):
	"""
	classdocs
	"""
	#
	#  Class methods
	#
	def __init__(self, proximity_structure, char_distribution, initial_acquaintance, 
	             interaction_likelihood_function,L,char_feedback=True,dynamic=True, degree_preference=None):
		"""
		Constructor
		"""
		#  Initialize instance variables
		self._proximity_structure = proximity_structure.copy()  #  Reference Proximity to which the model is relaxed to
		self._char_distribution = char_distribution.copy() #  Distribution of the individual characteristics
		self._char_feedback = char_feedback # flags if the character feedback is active or not
		adjacency=initial_acquaintance.adjacency
		self._acquaintance_network = Network(adjacency=adjacency, directed=False, silence_level=3)
		self._no_interactions=0.
		self._new_edges=0
		self._L=L
		self._char_weight=L.char_weight
		self._dynamic=dynamic


		if degree_preference is not None:
		    self._degree_preference = degree_preference.copy()
		else:
		    #  If no degree preference passed to __init__, use degree vector
		    #  of initial acquaintance network
		    self._degree_preference = initial_acquaintance.degree() 

		#  Set interaction_likelihood function
		self._interaction_likelihood_function = interaction_likelihood_function

		#  Set number of agents N
		self._N = self._acquaintance_network.N
		
		self._proximity_structure=proximity_structure
    
	#
	#  Service methods
	#

	def get_proximity_structure(self):
	    return self._proximity_structure

	def get_char_distribution(self):
	   return self._char_distribution

	def set_char_distribution(self,new_chd,dimension):
	   self._char_distribution[dimension,:]=new_chd

	def get_acquaintance_network(self):
	    return self._acquaintance_network

	def get_degree_preference(self):
	    return self._degree_preference
    
    #
    #  Model methods
    #
        
	def get_proximity_distances_characters(self,char_distribution,char_weight,proximity_structure):
		# returns a proximity matrix based on a distance between arbitrary distributed "individuals"
		distances = np.zeros((self._N,self._N))
		
		for i in xrange(self._N):
			#calculate the distances based on the first two chars that apply to all nodes uniformly
			distances[i,:] = char_weight[0]*np.abs(np.abs(char_distribution[0,:] - char_distribution[0,i])-1) + char_weight[1]*np.abs(np.abs(char_distribution[1,:] - char_distribution[1,i])-1) 
		return distances + char_weight[2]*proximity_structure

	def iterate(self, n_steps):
		"""
		Iterate the model for n_steps further time steps.
		"""

		#  Get current acquaintance network
		acquaintance = self.get_acquaintance_network()

		#  Get degree preference
		degree_preference = self.get_degree_preference()

		#  Get interaction_likelihood function
		interaction_likelihood_function = self._interaction_likelihood_function

		#  Iterate
		for i in xrange(n_steps):
			#  Calculate distance metric from acquaintance network
			distance_metric = self.get_distance_metric_matrix(acquaintance)

			#  Calculate interaction_likelihood matrix from distance metric matrix
			interaction_likelihood = self.get_interaction_likelihood_matrix(distance_metric,
			                                      interaction_likelihood_function)
			

			#  Obtain interaction_likelihood network by randomly drawing using interaction_likelihood
			#  probabilities
			interaction_likelihood_network = self.get_interaction_network(interaction_likelihood)


			# Interaction feedback on the individual characteristics
			# if dyn_meanfield==True, also for the dyn only case a probabalistic switch is applied
			if self._L.dyn_meanfield==True:
				char_distribution= self.update_char_feedback(self._char_distribution,interaction_likelihood_network)
			else:
				if self._char_feedback == True:
					char_distribution= self.update_char_feedback(self._char_distribution,interaction_likelihood_network)	
				else:
					char_distribution=self._char_distribution
					
			# Linear combination of the proximity derived via the characteristics and the proximity structure
			proximity =self.get_proximity_distances_characters(char_distribution,self._char_weight,self._proximity_structure)

			#  Update acquaintance network
			if self._dynamic==True:
				acquaintance = self.update_acquaintance_network(acquaintance, proximity,interaction_likelihood_network, degree_preference)
			
		
		#  Update instance variables to iterated values
		self._char_distribution = char_distribution
		self._acquaintance_network = acquaintance
		self._no_interactions = np.sum(interaction_likelihood_network)
		
	def get_eq_network(self):	
		"""
		Get the optimal target network for a given proximity by letting everybody interact with everybody else
		"""
		#  Get current acquaintance network
		acquaintance = self.get_acquaintance_network()

		#  Get degree preference
		degree_preference = self.get_degree_preference()
		
		# Get individual characteristics distribution
		char_distribution=self._char_distribution

		#  Obtain interaction_likelihood network by randomly drawing using interaction_likelihood
		#  probabilities
		interaction_likelihood_network = np.ones((self._N,self._N))

		# Linear combination of the proximity derived via the characteristics 
		proximity = self.get_proximity_distances_characters(char_distribution,self._char_weight,self._proximity_structure)
		
		#  Update acquaintance network
		acquaintance = self.update_acquaintance_network(acquaintance, proximity,interaction_likelihood_network, degree_preference)

		#  Update instance variables to iterated values
		self._char_distribution = char_distribution
		self._acquaintance_network = acquaintance
		self._no_interactions = np.sum(interaction_likelihood_network)		



        
	def get_distance_metric_matrix(self, acquaintance_network):
	    """
	    Simple variant: returns the geodesic graph distance (shortest path 
	    length) matrix
	    """
	    #  Return shortest path lengths matrix
	    return acquaintance_network.path_lengths()
    
	def get_interaction_likelihood_matrix(self, distance_metric_matrix, interaction_likelihood_function):
	    """
	    Returns the interaction_likelihood matrix given a distance metric on the 
	    acquaintance network using the given interaction_likelihood function.
	    """
	    return interaction_likelihood_function(distance_metric_matrix,self._L)
    
	def get_interaction_network(self, interaction_likelihood_matrix):
	    """
	    Returns the interaction_likelihood network's adjacency matrix given the interaction_likelihood matrix.
	    """
	    #  Draw uniformly distributed random numbers from the interval [0,1]
	    random_numbers = np.random.rand(self._N, self._N)
	    #  Symmetrize
	    random_numbers = (random_numbers + random_numbers.transpose()) / 2.
    
	    #  Return adjacency matrix of interaction_likelihood network
	    return (random_numbers <= interaction_likelihood_matrix).astype("int8")
    
	def update_char_feedback(self,char_distribution,interaction_likelihood_network):
		char_distribution_update=char_distribution
		smoking_matrix=self.get_acquaintance_network().node_attribute('smoker')
		sm_share=float(sum(smoking_matrix))/self._N
		#print 'sm_share',sm_share
		for i in xrange(char_distribution.shape[1]):
			nai=np.sum(interaction_likelihood_network[i,:])

			if nai ==0:
				char_distribution_update[:,i]= char_distribution[:,i]
			else:
				random_number=np.random.rand(1)
				# character feedback is based on social influence by the actual peers
				if self._char_feedback==True:
					if char_distribution[1,i]==0:
						prob_flip=.1*char_distribution[0,i]*np.sum(interaction_likelihood_network[i,:]*char_distribution[1,:])/nai
						char_distribution_update[1,i]=(random_number <= prob_flip).astype("int8")
					
					else: 	
						prob_flip=.1*(1-char_distribution[0,i])*(1-np.sum(interaction_likelihood_network[i,:]*char_distribution[1,:])/nai)
						char_distribution_update[1,i]=1-(random_number <= prob_flip).astype("int8")

				# mean field forcing for the dyn only case
				else:
					if self._L.dyn_meanfield:
						if char_distribution[1,i]==0:
							prob_flip=.1*char_distribution[0,i]*sm_share
							char_distribution_update[1,i]=(random_number <= prob_flip).astype("int8")
						
						else: 	
							prob_flip=.1*(1-char_distribution[0,i])*(1-sm_share)
							char_distribution_update[1,i]=1-(random_number <= prob_flip).astype("int8")



					
		# print np.sum(char_distribution_update[1,:]),np.sum(char_distribution_update[0,:])
		return char_distribution_update


	def update_acquaintance_network(self, acquaintance_network, 
	                                proximity_matrix, interaction_network,
	                                degree_preference):
	    """
	    Returns the updated acquantaince network's adjacency matrix.
	    """
		#  Get old acquaintance adjacency matrix
	    old_acquaintance_adjacency = acquaintance_network.adjacency
    
	    #  Get old degree k
	    k = acquaintance_network.degree()
    
	    #  Initialize new acquaintance adjacency matrix
	    acquaintance_adjacency = np.zeros(old_acquaintance_adjacency.shape, 
	                                      dtype="int8")
    
	    #  Initialize list of sorted indices
	    potential_acquaintance_indices = []
	    #  Loop over all agents 
	    for i in xrange(self._N):
	        #  Get node indices of acquaintances
	        acquaintance_indices = np.where(old_acquaintance_adjacency[i,:] == 1)[0]
        
	        if k[i] != 0:
	            #  Get node indices of interaction_likelihood_network
	            interaction_likelihood_indices = np.where(interaction_network[i,:] == 1)[0]
            
	            #  Combine both lists of indices and discard repeated entries
	            indices = list(np.unique(np.append(acquaintance_indices, 
	                                               interaction_likelihood_indices)))
            
	            #  Get proximity values of acquaintances and interaction_likelihood_network
	            similarities = proximity_matrix[i,indices]
            
	            indices = np.array(indices)
            
	            #  Get degree_preference[i] indices with largest similarities
	            sorted_indices = indices[similarities.argsort()][-degree_preference[i]:]  
            
	            #  Get k indices with largest similarities
	            potential_acquaintance_indices.append(sorted_indices)
               
	        else:
	            potential_acquaintance_indices.append([])
    
	    #  Keep only bidirectional potential new acquaintances
	    for i in xrange(self._N):
	        if k[i] != 0:
	            #  Initialize
	            filtered_indices = []
            
	            for j in potential_acquaintance_indices[i]:
	                if i in potential_acquaintance_indices[j]:
	                    filtered_indices.append(j)
            
	            acquaintance_adjacency[i,filtered_indices] = 1
        
	    #  Exclusively consider bidirectional acquaintance preferences
	    #acquaintance_adjacency = acquaintance_adjacency * acquaintance_adjacency.transpose()
    
	    #  Inclusively consider all bi- and unidirectional acquaintance preferences
	    #acquaintance_adjacency = (acquaintance_adjacency + acquaintance_adjacency.transpose() > 0).astype("int8")
	    
	    # return the no of changed changed contacts per turn (nodes can have edges with themselves, the diagonal is therefore removed)
	    self._new_edges= (np.sum(np.abs(old_acquaintance_adjacency-acquaintance_adjacency))-np.sum(np.diag(np.abs(old_acquaintance_adjacency-acquaintance_adjacency))))/2
		#print 'new edges:',np.sum((old_acquaintance_adjacency-acquaintance_adjacency))
    
	    #print "Symmetry:", ((acquaintance_adjacency - acquaintance_adjacency.transpose())**2).sum()
    
	    #  Update acquaintance Network object
	    acquaintance_network.adjacency=acquaintance_adjacency
            
	    #  Return new acquaintance network
	    return acquaintance_network
