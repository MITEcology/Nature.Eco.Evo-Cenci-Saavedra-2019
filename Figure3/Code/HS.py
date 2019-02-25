import heapq
import pylab as plt
from scipy import integrate
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import scipy
import numdifftools as nd
import sys, random
from matplotlib2tikz import save as tikz_save
from sklearn import preprocessing
################################################
### From A novel chaotic system with heart-shaped equilibrium and its circuital implementation
### Paham, Ciet-Thanh et al. Optik, 2016
################################################
f = open('hs.txt', 'w')

### Original parameters
a = 5.; b = 3.5;
#a = 2.; b = 8.; c = 2.95; k = 4.;
### Perturbed parameters


################################################
###### Initial conditions and time steps #######
x0 = .1; y0 = .1; z0 = .2;
T = 1000.;
dt = 0.01;
n_steps = T/dt;
t = np.linspace(0, T, n_steps)
X_f1 = np.array([x0, y0, z0])
######## Auxiliar Functions ####################
def dX_dt(X, t = 0):
	return(np.array([X[2],
			 -X[2]*(a*X[1] + b*X[1]**2 + X[0]*X[2]),
			X[0]**2 - abs(X[0])*X[1] + X[1]**2 -1]))
################################################
ts = integrate.odeint(dX_dt, X_f1, t)
X_f1 = ts[ts.shape[0]-1,:]
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
	####### Randomly perturbed
	K = random.randint(1,2)
	N = 2 - K
	sample = np.array([0] * N + [1] * K )
	np.random.shuffle(sample)
	sgm = random.uniform(5,10)

	a_ = a + np.random.normal(0, a/sgm)*sample[0];
	b_ = b + np.random.normal(0, b/sgm)*sample[1];

	def dX_dt_2(X, t = 0):
		return(np.array([X[2],
				 -X[2]*(a_*X[1] + b_*X[1]**2 + X[0]*X[2]),
				X[0]**2 - abs(X[0])*X[1] + X[1]**2 -1]))


	lower_percentile = np.percentile(trace,15)
	upper_percentile = np.percentile(trace,85)
	lower_indeces_div = [i for i,v in enumerate(trace) if v < lower_percentile]
	upper_indeces_div = [i for i,v in enumerate(trace) if v > upper_percentile]
	lower_indeces_div = np.asarray([i for i in lower_indeces_div if i >= 0])
	upper_indeces_div = np.asarray([i for i in upper_indeces_div if i >= 0])
	t_at_maximum_stability = random.choice(lower_indeces_div)*sampling_rate
	t_at_minimum_stability = random.choice(upper_indeces_div)*sampling_rate

	#print t_at_maximum_stability, t_at_minimum_stability, first_trace[idx_small_trace], first_trace[idx_large_trace]
	new_origin_at_min = np.squeeze(np.asarray(ts[t_at_maximum_stability,:]))
	new_origin_at_max = np.squeeze(np.asarray(ts[t_at_minimum_stability,:]))

	T = 1 #.5;
	dt = 0.01;
	n_steps = T/dt;
	t = np.linspace(0, T, n_steps)


	ts_at_min = integrate.odeint(dX_dt, new_origin_at_min, t)
	ts_at_max = integrate.odeint(dX_dt, new_origin_at_max, t)
	ts_perturbed_at_min = integrate.odeint(dX_dt_2, new_origin_at_min, t)
	ts_perturbed_at_max = integrate.odeint(dX_dt_2, new_origin_at_max, t)

	def distance(X,Y):
		return(np.sqrt(np.sum((X - Y)**2)))


	#### Normalize the time series
	min_max_scaler = preprocessing.MinMaxScaler()
	ts_at_min = min_max_scaler.fit_transform(ts_at_min)
	ts_at_max = min_max_scaler.fit_transform(ts_at_max)
	ts_perturbed_at_min = min_max_scaler.fit_transform(ts_perturbed_at_min)
	ts_perturbed_at_max = min_max_scaler.fit_transform(ts_perturbed_at_max)
	f.write('%i %lf %lf %lf\n' % (sample.sum(), distance(ts_at_min,ts_perturbed_at_min), distance(ts_at_max,ts_perturbed_at_max), np.log(distance(ts_at_min,ts_perturbed_at_min)/distance(ts_at_max,ts_perturbed_at_max))))
f.close()
