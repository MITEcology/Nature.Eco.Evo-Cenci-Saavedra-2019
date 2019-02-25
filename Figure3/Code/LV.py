import random, heapq
import pylab as plt
from scipy import integrate
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import scipy
import os, sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib2tikz import save as tikz_save
import numdifftools as nd
from sklearn import preprocessing
############################################################# Chaotic Lotka-Volterra ####################################
################################################
## open file
f = open('lv.txt', 'w')
########## Parameters adapted from "Chaos in low-dimensional ...."
r = np.array([1., 0.72, 1.53, 1.27])
A = np.matrix([[1.*r[0], 1.09*r[0], 1.52*r[0], 0.*r[0]],  [0.*r[1], 1.*r[1], 0.44*r[1],1.36*r[1]],  [2.33*r[2], 0.*r[2], 1.*r[2], 0.47*r[2]], [1.21*r[3], 0.51*r[3], 0.35*r[3], 1.*r[3]]])


#####################################################
#####################################################
################################################
################################################
###### Initial conditions and time steps #######
x0 = 0.2; y0 = 0.2; z0 = 0.3; k0 = 0.3;
####### At 1015 the model explode
T = 1000;
dt = 0.01;
n_steps = T/dt;
t = np.linspace(0, T, n_steps)
X_f1 = np.array([x0, y0, z0, k0])
################################################
################################################
def dX_dt(X, t = 0):
    dydt = np.array([X[s]*(r[s] - np.sum(np.dot(A,X)[0,s]))for s in range(0,len(X))])
    return(dydt)
################
################################################
printing = True
ini_cond = integrate.odeint(dX_dt, X_f1, t)
X_f1 = np.array([ini_cond[len(t)-1,0], ini_cond[len(t)-1,1], ini_cond[len(t)-1,2], ini_cond[len(t)-1,3]])

ts = integrate.odeint(dX_dt, X_f1, t)
trace = []
sampling_rate = 500
for i in range(0,ts.shape[0]):
	if i%500 == 0:
		f_jacob = nd.Jacobian(dX_dt)(np.squeeze(np.asarray(ts[i,:])))
		trace.append(np.trace(f_jacob))

#### min and max refer to the trace of the Jacobian
#### min trace == good structural stability
#### max trace == bad structural stability
number_of_iterations = int(sys.argv[1])
for it in range(number_of_iterations):
	K = random.randint(2,5)
	N = 16 - K
	sample = np.array([0] * N + [1] * K )
	np.random.shuffle(sample)
	sgm = random.uniform(10,15)

	A2 = np.matrix([[1.*r[0] + np.random.normal(0,A[0,0]/sgm)*sample[0], 1.09*r[0] + np.random.normal(0,A[0,1]/sgm)*sample[1], 1.52*r[0] + 
	np.random.normal(0,A[0,2]/sgm)*sample[2], 0.*r[0] + np.random.normal(0,A[0,3]/sgm)*sample[3]],
		       [0.*r[1]+ np.random.normal(0,A[1,0]/sgm)*sample[4], 1.*r[1]+ np.random.normal(0,A[1,1]/sgm)*sample[5], 0.44*r[1]+ 
	np.random.normal(0,A[1,2]/sgm)*sample[6],1.36*r[1] + np.random.normal(0,A[1,3]/sgm)*sample[7]],
		       [2.33*r[2]+np.random.normal(0,A[2,0]/sgm)*sample[8], 0.*r[2]+np.random.normal(0,A[2,1]/sgm)*sample[9], 
	1.*r[2]+np.random.normal(0,A[2,2]/sgm)*sample[10], 0.47*r[2] + np.random.normal(0,A[2,2]/sgm)*sample[11]],
		       [1.21*r[3]+np.random.normal(0,A[3,0]/sgm)*sample[12], 0.51*r[3]+np.random.normal(0,A[3,1]/sgm)*sample[13], 
	0.35*r[3]+np.random.normal(0,A[3,2]/sgm)*sample[14], 1.*r[3]+np.random.normal(0,A[3,3]/sgm)*sample[15]]])

	def dX_dt_2(X, t = 0):
	    dydt = np.array([X[s]*(r[s] - np.sum(np.dot(A2,X)[0,s]))for s in range(0,len(X))])
	    return(dydt)

	lower_percentile = np.percentile(trace,15)
	upper_percentile = np.percentile(trace,85)
	lower_indeces_div = [i for i,v in enumerate(trace) if v < lower_percentile]
	upper_indeces_div = [i for i,v in enumerate(trace) if v > upper_percentile]
	lower_indeces_div = np.asarray([i for i in lower_indeces_div if i >= 0])
	upper_indeces_div = np.asarray([i for i in upper_indeces_div if i >= 0])
	t_at_maximum_stability = random.choice(lower_indeces_div)*sampling_rate
	t_at_minimum_stability = random.choice(upper_indeces_div)*sampling_rate
	new_origin_at_min = np.squeeze(np.asarray(ts[t_at_maximum_stability,:]))
	new_origin_at_max = np.squeeze(np.asarray(ts[t_at_minimum_stability,:]))

	T = 20;
	dt = 0.01;
	n_steps = T/dt;
	t = np.linspace(0, T, n_steps)


	ts_at_min = integrate.odeint(dX_dt, new_origin_at_min, t)
	ts_at_max = integrate.odeint(dX_dt, new_origin_at_max, t)
	ts_perturbed_at_min = integrate.odeint(dX_dt_2, new_origin_at_min, t)
	ts_perturbed_at_max = integrate.odeint(dX_dt_2, new_origin_at_max, t)
	def distance(X,Y):
		return(np.sqrt(np.sum((X - Y)**2)))

	#### Normalize the time series to compare properly
	min_max_scaler = preprocessing.MinMaxScaler()
	ts_at_min = min_max_scaler.fit_transform(ts_at_min)
	ts_at_max = min_max_scaler.fit_transform(ts_at_max)
	ts_perturbed_at_min = min_max_scaler.fit_transform(ts_perturbed_at_min)
	ts_perturbed_at_max = min_max_scaler.fit_transform(ts_perturbed_at_max)
	f.write('%i %lf %lf %lf\n' % (sample.sum(), distance(ts_at_min,ts_perturbed_at_min), distance(ts_at_max,ts_perturbed_at_max), np.log(distance(ts_at_min,ts_perturbed_at_min)/distance(ts_at_max,ts_perturbed_at_max))))


f.close()
