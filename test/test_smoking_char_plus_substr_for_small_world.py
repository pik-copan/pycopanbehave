# Add repository path of the model core
import sys

sys.path.append('/Users/carls/Documents/PIK/ModelingSocialStability/Code/dynamic_society_model/src/model_core')
sys.path.append('/Users/carls/Documents/PIK/ModelingSocialStability/Code/src/model_core')
sys.path.append('/scratch/01/carls/social/emulation/model_core/')
sys.path.append('../../')
sys.path.append('../')
sys.path.append('../../../')
#  Import Network class from pygeonetwork
from pygeonetwork import Network

#  Import DynamicSocietyModel class from pydynamicsociety
from pysoc import DynamicSocietyModel

import igraph
#######################################
#  Number of nodes
N = 1000

#  Number of hysteresis iterations
n_transition = 1000

n_initial_eq_it=1
n_trans_eq_pre=0
n_trans_eq_post=50

n_iterations=n_transition+n_trans_eq_pre+n_trans_eq_post

# number of snapshots of the full nw
nw_snapshots=20
nw_save_steps=int(n_iterations/nw_snapshots)

#  Ensemble size (number of realizations of model time evolution)
n_ensemble = 20

#  Mean degree preference
mean_degree_pref = 10 # as it is stated in cf08
#  Degree preference standard deviation
std_degree_pref =3

# no of substructures in the Network
no_rings=30
links_per_node=2
# sets the mobility of the node in the substructure (relevant for creating the initial proximity matrix)
# three options: 'global': free mobility over the network, 'local': confined to substructure, 'ring': only between nearest neighbours in the initial underlying ring structure
mobility='global'

# additional to the initial structure, a number of nodes is connected randomly to ensure small world and high clustering at the same time
# the variable gives the no. of random connections relativ to the absolut no. of Nodes
# ws_random_reconnections=.1
# inter_substructure_random_connections=0.05

# strength of the transition into the less smoker state

alpha_trans=100*0.01

#probability of interaction
p_ai=.2
# offset that sets a basic interaction probability with all agents
interaction_offset=0.03


# sets the mobility with regard to the underlying proximity structure in positions
smoking_mobility=2
smoking_weight=.1*smoking_mobility
char_weight=(0,smoking_weight,(1-smoking_weight))



def interaction_likelihood_function(distance_metric_matrix):
	"""
	Returns the interaction_likelyhood matrix given a distance metric on the
	acquaintance network.
	The interaction likelyhood is treated exponential with a exp(-d/9)-relationship according 
	to the three degrees of influence 
	"""
	#  Example: Exponential relationship following fittings of by christakis and fowler plus plateau
	# return p_ai*((1.85-interaction_offset)* np.exp(-.6*distance_metric_matrix) + interaction_offset)
	
	# three degrees of influence plus an interaction offset to avoid disconnection
	return (p_ai-interaction_offset)*np.exp(-(distance_metric_matrix-1)/2.)+interaction_offset



def calc_cond_prob(smokers,nw_full,deg_sep_max):
	rcp=np.zeros(5)
	for i in xrange(deg_sep_max):
		deg_sep=i+1
		smoking_dep=[]
		for node in smokers:
			#print 'node',node
			distance_matrix=nw_full.get_path_lengths()	
			acquaintance_one=np.where(distance_matrix[node,:]==deg_sep)
			if acquaintance_one[0].size>0:
				smoking_dep.append(np.sum(nw_full.get_node_attribute('smoker')[acquaintance_one])/float(acquaintance_one[0].size)/(float(len(smokers))/float(N))-1)
		rcp[i]=np.mean(smoking_dep)
	return rcp

def distribution_function(yb,x):
	"""" modulate the parabolic distribution y= a(b-x)**2 +c fct. with the following boundary conditions
			1. x=0: y== 3 : follows from the initial distribution function
			2. Int[0,1]==1: normalization criteria
			3. x=1: y==yb, where yb is the tunable parameter, which reflects the probability for an agent to have smoking disposition 1 
			yb0==3
	"""
	# print 'yb',yb
	#yb=3
	b=(1.+yb/3)/(1+yb)
	a=2./(b-1./3)
	c=3.-a*b**2

	return   a*(b-x)**2 +c

def rejection_sampling(N, distribution_function_ini,yb):
	#  Creates random sampling of an arbitraty distribution using the rejection sampling method for the continous characteristic (eg. smoking affinity) and a second binary characteristic based on this
	result = np.zeros((2,N))
	i = 0
	while i < N:
		random_numbers = np.random.rand(2)
		if random_numbers[0] < distribution_function_ini(yb,random_numbers[1]):
			result[0,i] = random_numbers[1]
			# if the smoking preference is greater than a certain value, the binary characteristic smoker/non-smoker is assigned
			result[1,i]=(random_numbers[1] > .5).astype('int8')
			i+=1
	return result

def proximity_str_ring(substructure1):
	small_world_distance_matrix=np.asarray(substructure1.shortest_paths())
	ring_length=small_world_distance_matrix.shape[0]
	distance=1-.1*(small_world_distance_matrix-1)
	distance[np.where(distance<=0)]=0
	for k in xrange(distance.shape[0]):
		distance[k,k]=0

	for i in xrange(ring_length):
		# find the top edges for the individual agent i
		dist_i=distance[i,:]
		top_edges=dist_i.argsort()[-2*links_per_node:]
		for l in top_edges:
			if np.random.random() <= ws_random_reconnections:
				m,n=ring_length*np.random.random(2)
				m,n=int(m),int(n)
				while n==m:
					n=int(ring_length*np.random.random())
				old_edge=dist_i[l].copy()
				new_edge=distance[m,n].copy()			
				distance[i,l]=distance[l,i]=new_edge
				distance[m,n]=distance[n,m]=old_edge
				#print 'random reconnected',i,l,'with',m,n	
	return distance

def generate_initial_distance_sm(mobility,inter_substructure_random_connections):
	#inter_substructure_random_connections=.5
	
	substr_adj_list=[]
	ringlength_standard=N/no_rings#np.asarray(np.linspace(0,N,no_rings+1),dtype='int')
	full_adjacency=np.zeros((N,N))
	proximity_small_world=np.zeros((N,N))
	for ring in xrange(no_rings):
		ringlength=ringlength_standard
		if ring == (no_rings-1):
			#print 'last ring'
			#account for the last nodes missing
			ringlength=N-ringlength_standard*(no_rings-1)-1
		substructure1=igraph.GraphBase.Ring(ringlength)
		
		substructure1=igraph.GraphBase.Lattice([ringlength], nei=links_per_node, directed=False, mutual=True, circular=True)
		if mobility !='ring':
			substructure1.rewire(int(ws_random_reconnections*ringlength*links_per_node))
		
		elif mobility =='ring':
			proximity_small_world[ring*ringlength_standard:ring*ringlength_standard+ringlength,ring*ringlength_standard:ring*ringlength_standard+ringlength]=proximity_str_ring(substructure1)
		
		substr_adj_list.append(np.asarray(substructure1.get_adjacency()))
		full_adjacency[ring*ringlength_standard:ring*ringlength_standard+ringlength,ring*ringlength_standard:ring*ringlength_standard+ringlength]=np.asarray(substructure1.get_adjacency())
	
	adj_list=[]
	for k in xrange(N):
		adj_list.append(full_adjacency[k,:])
	full_structure=igraph.GraphBase.Adjacency(adj_list, mode='undirected')
	if mobility =='global':
		full_structure.rewire(int(inter_substructure_random_connections*N*links_per_node))
		print 'apl initial NW',full_structure.average_path_length()
		small_world_distance_matrix=np.asarray(full_structure.shortest_paths())
		proximity_small_world=1-.1*(small_world_distance_matrix-1)
		proximity_small_world[np.where(proximity_small_world<=0)]=0
	elif mobility =='local':
		small_world_distance_matrix=np.asarray(full_structure.shortest_paths())
		proximity_small_world=1-.1*(small_world_distance_matrix-1)
		proximity_small_world[np.where(proximity_small_world<=0)]=0
	if 	mobility !='global':
		for k in xrange(int(inter_substructure_random_connections*N*links_per_node)):
			rn1=np.random.randint(N)
			rn2=np.random.randint(N)
			proximity_small_world[rn1,rn2]=proximity_small_world[rn2,rn1]=1
			
	for k in xrange(proximity_small_world.shape[0]):
		proximity_small_world[k,k]=0
		#imshow(proximity_small_world)
		#title('proxim structure, int rc %s, apl %s'%(inter_substructure_random_connections,full_structure.average_path_length()))
		
	return proximity_small_world
	# return full_structure.average_path_length(), full_structure.transitivity_avglocal_undirected()


apl=[]
clustering=[]
cond_prob=[]
inter_substructure_random_connections_range=np.linspace(0.01,0.2,20)
ws_random_reconnections=.01
for inter_substructure_random_connections in inter_substructure_random_connections_range:

	char_distribution_initial= rejection_sampling(N, distribution_function,3)

	#  Degree preference distribution
	degree_preference = np.random.normal(mean_degree_pref, std_degree_pref, N).astype("int8")
	acquaintance_network = Network.ErdosRenyi(n_nodes=N, link_probability=mean_degree_pref / float(N - 1))
	proximity_structure=generate_initial_distance_sm(mobility,inter_substructure_random_connections)
	char_feedback=True
	dynamic=True
	model_initial=DynamicSocietyModel (proximity_structure,char_distribution_initial, acquaintance_network,
	                 interaction_likelihood_function,char_weight,char_feedback,dynamic,degree_preference)
	model_initial.get_eq_network()
	model_initial.iterate(10)
	model_initial.get_acquaintance_network().set_node_attribute('smoker',model_initial._char_distribution[1,:].copy())
	apl.append(model_initial.get_acquaintance_network().get_average_path_length())
	clustering.append(model_initial.get_acquaintance_network().get_global_clustering())
	cond_prob.append(calc_cond_prob(np.where(model_initial.get_char_distribution()[1,:]==1)[0],model_initial.get_acquaintance_network(),1))
	#print substructure1.average_path_length()

figure()
plt.subplot(2,1,1)
plt.plot(inter_substructure_random_connections_range,apl,marker='x',label='apl')
plt.legend(loc='best')
plt.ylabel('apl')
plt.xlabel('rewiring prob')
plt.twinx()
plt.plot(inter_substructure_random_connections_range,clustering,marker='o',label='clustering',color='green')
plt.legend(loc='best')
plt.ylabel('clustering')
plt.title('small-world properties as a function of the inner-substructure rewiring \n ws_random_reconnections %s sm mobility: 2'%ws_random_reconnections)
plt.xlabel('rewiring prob')

plt.subplot(2,1,2)
plt.plot(inter_substructure_random_connections_range,np.asarray(cond_prob)[:,0],marker='o',label='cond_prob',color='green')
plt.legend(loc='best')
plt.xlabel('rewiring prob')
plt.ylabel('initial cond probability')
plt.savefig('test_sm_char_plus_substr_wsr_%s.pdf'%ws_random_reconnections)
#plt.errorbar(ws_random_reconnections_range,clustering/clustering.max(),label='C(p)/C(0)',marker='o')

#imshow(np.asarray(substructure1.get_adjacency()))