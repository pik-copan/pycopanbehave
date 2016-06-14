#!/usr/bin/python

"""" 		
	this script derives and plots the mean field number of smokers over a range  of different target disposition distributions
"""



import pylab as plt
import numpy as np

def final_distribution_function(yb,x):
	"""" modulate the parabolic distribution y= a(b-x)**2 +c fct. with the following boundary conditions
			1. x=0: y== 3 : follows from the initial distribution function
			2. Int[0,1]==1: normalization criteria
			3. x=1: y==yb, where yb is the tunable parameter, which reflects the probability for an agent to have smoking disposition 1 
			yb0==3
	"""
	b=(1.+yb/3)/(1+yb)
	a=2./(b-1./3)
	c=3.-a*b**2
	#print 'yb in function', yb
	return   a*(b-x)**2 +c 	


def rejection_sampling(N, distribution_function_ini):
	"""
	Creates random sampling of an arbitraty distribution using the rejection sampling method 
	for the continous characteristic (eg. smoking affinity) and a second binary characteristic based on this.
	
	Returns a 2xN array.
	"""
	result = np.zeros((2,N))
	i = 0
	while i < N:
		random_numbers = np.random.rand(2)
		if random_numbers[0] < distribution_function_ini(yb,random_numbers[1]):
			result[0,i] = random_numbers[1]
			# if the smoking preference is greater than a certain value, the binary characteristic smoker/non-smoker is assigned
			if random_numbers[1] > .5:
							result[1,i]=1
			result[1,i]=(random_numbers[1] > .5).astype('int8')
			i+=1
	return result


	
N=1000
ens_size=10
range_size=20
range_yb=np.linspace(0,3,range_size)
iteration_steps=20
output_matrix=np.zeros((range_size,iteration_steps,ens_size))


for l in xrange(ens_size):
	i=0
	for yb in range_yb:
		# the target distribution is created
		
		rej_semp=rejection_sampling(N, final_distribution_function)
		
		# the system is initialized with the individual only expectation value (threshold .5 in the disposition distribution)
		output_matrix[i,0,l]=sum(rej_semp[1,:])
		
		for k in xrange(1,iteration_steps):
			sm_t_min_1=float(output_matrix[i,k-1,l])
			
			try:
				output_matrix[i,k,l]=np.where(rej_semp[0,:]/(1-rej_semp[0,:])/(N/sm_t_min_1-1)>1)[0].shape[0]
			
			except:'dev by 0'
				#print 'dev by 0'
		i+=1
		
mean_field_std=output_matrix.std(axis=2)[:,-1]
mean_field_mean=output_matrix.mean(axis=2)[:,-1]

ind_expect_std=output_matrix.std(axis=2)[:,0]
ind_expect_mean=output_matrix.mean(axis=2)[:,0]

plt.figure()
plt.errorbar(range_yb,ind_expect_mean,yerr=ind_expect_std,label='Ind expect',color='green' )
plt.errorbar(range_yb,mean_field_mean,yerr=mean_field_std,label='Mean Field',color='black' )
plt.title('No of smokers vs. disp distr')
plt.xlabel('y(x==1)')
plt.ylabel('smokers')
plt.legend(loc='best')
plt.savefig('no_smokers_mean_field.pdf')
