import igraph
N=1000
ringlength=33

ws_random_reconnections=.1
apl=[]
clustering=[]
ws_random_reconnections_range=np.logspace(-3,-.5,100)
for ws_random_reconnections in ws_random_reconnections_range:
	# ws_random_reconnections=.01
	substructure1=igraph.GraphBase.Ring(ringlength)
	# substructure1.rewire(int(ws_random_reconnections*ringlength))
	# #apl.append(substructure1.average_path_length())
	# print substructure1.average_path_length()
	#ws_random_reconnections=.01
	links_per_node=2
	substructure1=igraph.GraphBase.Lattice([33], nei=links_per_node, directed=False, mutual=True, circular=True)
	substructure1.rewire(int(ws_random_reconnections*ringlength*links_per_node))
	
	apl.append(substructure1.average_path_length())
	clustering.append(np.asarray(substructure1.transitivity_local_undirected()).mean())
	#print substructure1.average_path_length()

figure()
plt.plot(ws_random_reconnections_range,apl,marker='x',label='apl')
plt.legend(loc='best')
plt.twinx()
plt.plot(ws_random_reconnections_range,clustering,marker='o',label='clustering',color='green')
plt.legend(loc='best')
plt.xlabel('rewiring prob')

plt.savefig('apl_vs_recon_substructure.pdf')
#plt.errorbar(ws_random_reconnections_range,clustering/clustering.max(),label='C(p)/C(0)',marker='o')

#imshow(np.asarray(substructure1.get_adjacency()))