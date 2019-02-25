import heapq
import pylab as plt
from scipy import integrate,stats
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import scipy
import numdifftools as nd
import sys, random
from matplotlib2tikz import save as tikz_save
from sklearn import preprocessing

################################################
## open file
f = open('CR.txt', 'w')
##

nu1 = 0.1; nu2 = 0.07;
C1 = 0.5; C2 = 0.5;
lambda1 = 3.2; lambda2 =2.9;
mu1 = 0.15; mu2 = 0.15;
kappa1 = 2.5; kappa2 = 2.;
Rstar = 0.3; k = 1.2;
################################################
###### Initial conditions and time steps #######
### These initial conditions are on the attractor
p1_0 = 0.006884225; p2_0 = 0.087265965; c1_0 = 0.002226393; c2_0 = 1.815199890; r_0 = 0.562017616;
T = 5000.;
dt = 0.01;
n_steps = T/dt;
t = np.linspace(0, T, n_steps)
X_f1 = np.array([p1_0, p2_0, c1_0, c2_0, r_0])
######## Auxiliar Functions ####################
def Uptake(var_x, L, KI):
    return(L*var_x/(KI + var_x))
################################################
def dX_dt(X, t = 0):
    dydt = np.array([nu1*Uptake(X[2], lambda1, C1)*X[0] - nu1*X[0],
                     nu2*Uptake(X[3], lambda2, C2)*X[1] - nu2*X[1],
                     mu1*Uptake(X[4], kappa1, Rstar)*X[2] - mu1*X[2] - nu1*Uptake(X[2], lambda1, C1)*X[0],
                     mu2*Uptake(X[4], kappa2, Rstar)*X[3] - mu2*X[3] - nu2*Uptake(X[3], lambda2, C2)*X[1],
                     X[4]*(1 - X[4]/k) -  mu1*Uptake(X[4], kappa1, Rstar)*X[2] - mu2*Uptake(X[4], kappa2, Rstar)*X[3]])
    return(dydt)


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
	K = random.randint(2,8)
	N = 12 - K
	sample = np.array([0] * N + [1] * K )
	np.random.shuffle(sample)
	sgm = random.uniform(5,10)
	nu1_ = 0.1 + np.random.normal(0,nu1/sgm)*sample[0]; nu2_ = 0.07 + np.random.normal(0,nu2/sgm)*sample[1];
	C1_ = 0.5 + np.random.normal(0,C1/sgm)*sample[2]; C2_ = 0.5 + np.random.normal(0,C2/sgm)*sample[3];
	lambda1_ = 3.2 + np.random.normal(0,lambda1/sgm)*sample[4]; lambda2_ =2.9 + np.random.normal(0,lambda2/sgm)*sample[5];
	mu1_ = 0.15 + np.random.normal(0,mu1/sgm)*sample[6]; mu2_ = 0.15 + np.random.normal(0,mu2/sgm)*sample[7];
	kappa1_ = 2.5 + np.random.normal(0,kappa1/sgm)*sample[8]; kappa2_ = 2. + np.random.normal(0,kappa2/sgm)*sample[9];
	Rstar_ = 0.3 + np.random.normal(0,Rstar/sgm)*sample[10]; k_ = 1.2 + np.random.normal(0,k/sgm)*sample[11];
	def dX_dt_2(X, t = 0):
	    dydt = np.array([nu1_*Uptake(X[2], lambda1_, C1_)*X[0] - nu1_*X[0],
	                     nu2_*Uptake(X[3], lambda2_, C2_)*X[1] - nu2_*X[1],
	                     mu1_*Uptake(X[4], kappa1_, Rstar_)*X[2] - mu1_*X[2] - nu1_*Uptake(X[2], lambda1_, C1_)*X[0],
	                     mu2_*Uptake(X[4], kappa2_, Rstar_)*X[3] - mu2_*X[3] - nu2_*Uptake(X[3], lambda2_, C2_)*X[1],
	                     X[4]*(1 - X[4]/k_) -  mu1_*Uptake(X[4], kappa1_, Rstar)*X[2] - mu2_*Uptake(X[4], kappa2_, Rstar_)*X[3]])
	    return(dydt)

	lower_percentile = np.percentile(trace,10)
	upper_percentile = np.percentile(trace,90)
	lower_indeces_div = [i for i,v in enumerate(trace) if v < lower_percentile]
	upper_indeces_div = [i for i,v in enumerate(trace) if v > upper_percentile]
	lower_indeces_div = np.asarray([i for i in lower_indeces_div if i >= 0])
	upper_indeces_div = np.asarray([i for i in upper_indeces_div if i >= 0])
	t_at_maximum_stability = random.choice(lower_indeces_div)*sampling_rate
	t_at_minimum_stability = random.choice(upper_indeces_div)*sampling_rate

	new_origin_at_min = np.squeeze(np.asarray(ts[t_at_maximum_stability,:]))
	new_origin_at_max = np.squeeze(np.asarray(ts[t_at_minimum_stability,:]))

	T = 40;
	dt = 0.01;
	n_steps = T/dt;
	t = np.linspace(0, T, n_steps)


	ts_at_min = integrate.odeint(dX_dt, new_origin_at_min, t)
	ts_at_max = integrate.odeint(dX_dt, new_origin_at_max, t)
	ts_perturbed_at_min = integrate.odeint(dX_dt_2, new_origin_at_min, t)
	ts_perturbed_at_max = integrate.odeint(dX_dt_2, new_origin_at_max, t)
	def distance(X,Y):
		return(np.sqrt(np.sum((X - Y)**2)))
	#### Standardize time series
	min_max_scaler = preprocessing.MinMaxScaler()
	ts_at_min = min_max_scaler.fit_transform(ts_at_min)
	ts_at_max = min_max_scaler.fit_transform(ts_at_max)
	ts_perturbed_at_min = min_max_scaler.fit_transform(ts_perturbed_at_min)
	ts_perturbed_at_max = min_max_scaler.fit_transform(ts_perturbed_at_max)
	#################		maximum stability 			minimum stability					Ratio
	f.write('%i %lf %lf %lf\n' % (sample.sum(), distance(ts_at_min,ts_perturbed_at_min), distance(ts_at_max,ts_perturbed_at_max), np.log(distance(ts_at_min,ts_perturbed_at_min)/distance(ts_at_max,ts_perturbed_at_max))))





f.close()
