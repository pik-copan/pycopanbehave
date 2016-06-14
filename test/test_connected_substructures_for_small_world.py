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
#  Number of nodes
N = 1000
# no of substructures in the Network
no_rings=30
links_per_node=2
# sets the mobility of the node in the substructure (relevant for creating the initial proximity matrix)
# three options: 'global': free mobility over the network, 'local': confined to substructure, 'ring': only between nearest neighbours in the initial underlying ring structure
mobility='global'


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
		
	# return proximity_small_world
	return full_structure.average_path_length(), full_structure.transitivity_avglocal_undirected()

ringlength=33

ws_random_reconnections=.1
apl=[]
clustering=[]
inter_substructure_random_connections_range=np.logspace(-3,0,20)
for inter_substructure_random_connections in inter_substructure_random_connections_range:
	apl1,cl1=generate_initial_distance_sm('global',inter_substructure_random_connections)
	apl.append(apl1)
	clustering.append(cl1)
	#print substructure1.average_path_length()

figure()
plt.plot(inter_substructure_random_connections_range,apl,marker='x',label='apl')
plt.legend(loc='best')
plt.ylabel('apl')
plt.xlabel('rewiring prob')
plt.twinx()
plt.plot(inter_substructure_random_connections_range,clustering,marker='o',label='clustering',color='green')
plt.legend(loc='best')
plt.ylabel('clustering')
plt.title('small-world properties as a function of the inner-substructure rewiring \n ws_random_reconnections %s'%ws_random_reconnections)
plt.xlabel('rewiring prob')

plt.savefig('apl_vs_recon_fullstructure_wsr_%s.pdf'%ws_random_reconnections)
#plt.errorbar(ws_random_reconnections_range,clustering/clustering.max(),label='C(p)/C(0)',marker='o')

#imshow(np.asarray(substructure1.get_adjacency()))