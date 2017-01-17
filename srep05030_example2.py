import numpy as np
import math

from srep05030_network_reconstruction import network_reconstruction


def example2():

	#----------------------------------------
	#The connection matrix as shown in Figure 1a of the paper
	#----------------------------------------
	A0 = [[0.0,0.7,0.32,0.0,0.0,0.04],[0.7,0.0,0.9,0.0,0.0,0.0],[0.845,0.91,0.75,0.11,0.0,0.0],[0.13,0.6,0.68,0.0,0.22,0.0],[0.69,0.0,0.0,0.0,0.0,0.0],[0.0,0.27,0.0,0.95,0.375,0.0]]

	#----------------------------------------
	#This is the original time series data Figure 2, note that the initial conditions are randomly generated (from a normal distribution) so some of the values may vary
	#----------------------------------------
	t = [0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0]
	x = [   [0.5688,0.4694,0.0119,0.3371,0.1622,0.7943],
		 	[0.4700,0.3889,0.0182,0.5690,0.2533,0.9345],
	 		[0.3890,0.3202,0.0205,0.7603,0.3307,1.0324],
			[0.3229,0.2629,0.0240,0.9173,0.3951,1.0819],
			[0.2689,0.2155,0.0297,1.0459,0.4483,1.0910],
			[0.2245,0.1766,0.0368,1.1508,0.4920,1.0744],
			[0.1880,0.1446,0.0444,1.2363,0.5278,1.0443],
			[0.1578,0.1184,0.0517,1.3058,0.5572,1.0085],
			[0.1327,0.0970,0.0585,1.3623,0.5813,0.9715],
			[0.1118,0.0794,0.0646,1.4081,0.6010,0.9358],
			[0.0944,0.0650,0.0699,1.4453,0.6171,0.9027],
			[0.0799,0.0532,0.0746,1.4753,0.6303,0.8726],
			[0.0677,0.0436,0.0786,1.4996,0.6411,0.8456],
			[0.0575,0.0357,0.0821,1.5192,0.6500,0.8219],
			[0.0490,0.0292,0.0850,1.5350,0.6572,0.8010],
			[0.0418,0.0239,0.0875,1.5478,0.6632,0.7829],
			[0.0358,0.0196,0.0896,1.5581,0.6680,0.7672],
			[0.0308,0.0160,0.0913,1.5664,0.6720,0.7536],
			[0.0266,0.0131,0.0928,1.5730,0.6753,0.7420],
			[0.0231,0.0107,0.0941,1.5784,0.6780,0.7320],
			[0.0201,0.0088,0.0951,1.5827,0.6801,0.7235] ]

	
	#Plot the time series by using
	
	import matplotlib.pylab as plt
	plt.plot(np.matrix(t).T,np.matrix(x))
	plt.xlabel('Time')
	plt.ylabel('Gene regulation')
	plt.savefig('example2_time_series.png')
	

	#----------------------------------------
	#the function for example 2 as shown in Eq.1 of the paper
	#----------------------------------------
	def f(x):
		return - x
	def h(x,i,k):
		if k < 3:
			return x**5 / (1.0+x**5)
		else:
			return 1.0 / (1.0+x**5)

	reconstructed_A = network_reconstruction(x,A0,f,h,R=10)

		
	#Visualize the original adjacency matrix and the predicted one and compare you can use
	
	import matplotlib.pylab as plt
	fig, axes = plt.subplots(nrows=1, ncols=3)


	plt.set_cmap('bwr')
	vmin,vmax = -1,1
	im = axes.flat[0].imshow(A0,interpolation='none', vmin=vmin, vmax=vmax)
	axes.flat[0].set_title('Original A')
	plt.set_cmap('bwr')
	im = axes.flat[1].imshow(reconstructed_A,interpolation='none', vmin=vmin, vmax=vmax)
	axes.flat[1].set_title('Reconstructed A')
	plt.set_cmap('bwr')
	im = axes.flat[2].imshow(np.abs(reconstructed_A-A0),interpolation='none', vmin=vmin, vmax=vmax)
	axes.flat[2].set_title('Difference')

	fig.subplots_adjust(right=0.8)
	cbar_ax = fig.add_axes([0.85, 0.35, 0.02, 0.3])
	fig.colorbar(im, cax=cbar_ax)
	plt.savefig('example2_adjacency_matrix.png')
	

if __name__ == "__main__":
	example2()


