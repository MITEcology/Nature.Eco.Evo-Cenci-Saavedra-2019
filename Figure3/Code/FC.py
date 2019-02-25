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
## open file
f = open('FC.txt', 'w')
### Original parameters
r = 4.3; k = 50; a1 = b1 = 0.1; s = 1.; h = 0.05; a2 = b2 = 0.1; l = 1; n = 0.03

################################################
###### Initial conditions and time steps #######
#p1_0 = 0.7; p2_0 = 0.8; c1_0 = 0.5; c2_0 = 0.8; r_0 = 1.;
p1_0 = 0.006884225; p2_0 = 0.087265965; c1_0 = 0.002226393; c2_0 = 1.815199890; r_0 = 0.562017616;
T = 1000.;
dt = 0.01;
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
	K = random.randint(2,5)
	N = 9 - K
	sample = np.array([0] * N + [1] * K )
	np.random.shuffle(sample)
	sgm = random.uniform(6,10)

	### Perturbed parameters
	r_ = 4.3 # + np.random.normal(0,0.01)*sample[9];
	k_ = 50 + np.random.normal(0,k/sgm)*sample[0]; a1_ = 0.1 + np.random.normal(0,a1/sgm)*sample[1]
	b1_ = 0.1 + np.random.normal(0,b1/sgm)*sample[2]; s_ = 1.+ np.random.normal(0,s/sgm)*sample[3]; h_ = 0.05+ np.random.normal(0,h/sgm)*sample[4]; 
	a2_ = 0.1+ np.random.normal(0,a2/sgm)*sample[5]
	b2_ = 0.1+ np.random.normal(0,b2/sgm)*sample[6];
	l_ = 1+ np.random.normal(0,l/sgm)*sample[7]; n_ = 0.03+ np.random.normal(0,n/sgm)*sample[8]

	def dX_dt_2(X, t = 0):
		return(np.array([r_*X[0]*(1-X[0]/k_) - a1_*X[0]*X[1]/(1+b1_*X[0]),
				-s_*X[1] + h_*X[0]*X[1] - a2_*X[1]*X[2]/(1+b2_*X[1]),
				-l_*X[2] + n_*X[1]*X[2]]))

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

	T = 5;
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
