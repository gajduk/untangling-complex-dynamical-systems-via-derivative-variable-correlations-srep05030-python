import numpy as np
import math
import random

def network_reconstruction(x,A0,f,h,n_tries=1000,tau=0.2,R=1):
	'''
	:param x - time series data, a matrix with N columns and L rows, where N is the number of proteins and L is the number of time points
	:param A0 - adjacency matrix N rows and N columns
	:param f - see Eq.1 in the paper - accepts a single double number, and returns a single double number e.g. f(x) return -x
	:param h - see Eq.1 in the paper - accepts a single double number, and returns a single double number e.g. h(x) return math.tanh(x)
	:param tau - time between two measurements, all time points should be equally spaced
	:param R - not sure what it is
	return - the reconstructed adjacency matrix
	'''

	def g_generator(x_min,x_max,n_harmonics=10):
		pi = np.pi
		def coef_generator():
			return math.pow(10, random.random() * 2.00432137) - 1.0
		a,b = np.matrix([coef_generator() for _ in range(n_harmonics)]),np.matrix([coef_generator() for _ in range(n_harmonics)])
		k = np.matrix([i + 1 for i in range(n_harmonics)])
		def g(x):
			x_norm = k*pi*(x-x_min)/(x_max-x_min)
			return  np.sum(np.multiply(a,np.sin(x_norm)) + np.multiply(b,np.cos(x_norm)))
		return g

	#L is the length of the time series, 
	L, N = np.matrix(x).shape
	x = np.matrix(x).T.tolist()
	xx,dx,f_F,h_F = np.zeros((N,L-1)),np.zeros((N,L-1)),np.zeros((N,L-1)),np.zeros((N,L-1))

	for i in range(N):
		for k in range(L-1):
			xx[i][k] = ( x[i][k] + x[i][k+1] ) / 2.0
			dx[i][k] = ( x[i][k+1] - x[i][k] )
			f_F[i][k] = f(xx[i][k])
			h_F[i][k] = h(xx[i][k])
	min_delta = 100000
	reconstructed_A = -1
	deltas = []
	x_min,x_max = np.amin(x),np.amax(x)
	for iter_ in range(n_tries):
		g = g_generator(x_min,x_max)
		g_F = np.zeros((N,L-1))

		for i in range(N):
			for k in range(L-1):
				g_F[i][k] = g(xx[i][k])

		B, C, E = np.zeros((N,N)),np.zeros((N,N)),np.zeros((N,N))
		for i in range(N):
			for j in range(N):
				B[i][j] = np.sum(np.multiply(g_F[:][i],dx[:][j])) / ( tau * (L-1) * R)
				C[i][j] = np.sum(np.multiply(g_F[:][i],f_F[:][j]))  / ( 1.0 * (L-1) * R)
				E[i][j] = np.sum(np.multiply(g_F[:][i],h_F[:][j])) / ( 1.0 * (L-1) * R)

		A = np.dot(np.linalg.inv(E),(B - C))

		delta_A = 0.0
		normalization = 0
		for i in range(N):
			for j in range(N):
				delta_A += (A[i][j]-A0[i][j])**2
				normalization += A0[i][j]**2
		delta_A = math.sqrt(delta_A/normalization)
		deltas.append(delta_A)
		if delta_A < min_delta:
			min_delta = delta_A
			reconstructed_A = A
			
	return reconstructed_A