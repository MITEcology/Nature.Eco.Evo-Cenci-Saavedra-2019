import pylab as plt
from scipy import integrate
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import scipy
import numdifftools as nd
import sys, random
### Original parameters
r = 4.3; k = 50; a1 = b1 = 0.1; s = 1.; h = 0.05; a2 = b2 = 0.1; l = 1; n = 0.03
###### Initial conditions and time steps #######
p1_0 = 0.006884225; p2_0 = 0.087265965; c1_0 = 0.002226393;
T = 500.;
dt = 0.001;
n_steps = T/dt;
t = np.linspace(0, T, n_steps)
X_f1 = np.array([p1_0, p2_0, c1_0])
######## Auxiliar Functions ####################
def dX_dt(X, t = 0):
	return(np.array([r*X[0]*(1-X[0]/k) - a1*X[0]*X[1]/(1+b1*X[0]),
			-s*X[1] + h*X[0]*X[1] - a2*X[1]*X[2]/(1+b2*X[1]),
			-l*X[2] + n*X[1]*X[2]]))
################################################
ts = integrate.odeint(dX_dt, X_f1, t)
X_f1 = ts[ts.shape[0]-1,:]
ts = integrate.odeint(dX_dt, X_f1, t)

g = open('../inputFiles/deterministic_chaos_fc.txt', 'w')
jacobian_matrix = open('../inputFiles/jacobian_chaos_fc.txt', 'w')

num_species = ts.shape[1]
for i in range(0,ts.shape[0]):
	if i%200 == 0:
		f_jacob = nd.Jacobian(dX_dt)(np.squeeze(np.asarray(ts[i,:])))
		g.write('%f %f %f\n' % (ts[i, 0], ts[i, 1], ts[i, 2]))
		for u in range(0,num_species):
			for z in range(0,num_species):
			    	jacobian_matrix.write('%lf ' % (f_jacob[u,z]))
		jacobian_matrix.write('\n')
g.close()
jacobian_matrix.close()
