import pylab as plt
from scipy import integrate
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import scipy
import numdifftools as nd

#### Code to generate the time series of species abundance for the 5 species consumer-resources model
#### Here we also compute the Jacobian coefficients analytically so that we can then check
#### If the VCR is correctly inferred from the regularized S-map

#### run as python CR.py

###### Define the parameters
nu1 = 0.1; nu2 = 0.07;
C1 = 0.5; C2 = 0.5;
lambda1 = 3.2; lambda2 =2.9;
mu1 = 0.15; mu2 = 0.15;
kappa1 = 2.5; kappa2 = 2.;
Rstar = 0.3; k = 1.2;
###### Initial conditions and time steps #######
p1_0 = 0.006884225; p2_0 = 0.087265965; c1_0 = 0.002226393; c2_0 = 1.815199890; r_0 = 0.562017616;
T = 5000.;
dt = 0.01;
n_steps = T/dt;
t = np.linspace(0, T, n_steps)
X_f1 = np.array([p1_0, p2_0, c1_0, c2_0, r_0])
######## Model
def Uptake(var_x, L, KI):
    return(L*var_x/(KI + var_x))
def dX_dt(X, t = 0):
    dydt = np.array([nu1*Uptake(X[2], lambda1, C1)*X[0] - nu1*X[0],
                     nu2*Uptake(X[3], lambda2, C2)*X[1] - nu2*X[1],
                     mu1*Uptake(X[4], kappa1, Rstar)*X[2] - mu1*X[2] - nu1*Uptake(X[2], lambda1, C1)*X[0],
                     mu2*Uptake(X[4], kappa2, Rstar)*X[3] - mu2*X[3] - nu2*Uptake(X[3], lambda2, C2)*X[1],
                     X[4]*(1 - X[4]/k) -  mu1*Uptake(X[4], kappa1, Rstar)*X[2] - mu2*Uptake(X[4], kappa2, Rstar)*X[3]])
    return(dydt)
################
ts = integrate.odeint(dX_dt, X_f1, t)
first_trace = []
g = open('../inputFiles/deterministic_chaos_cr.txt', 'w')
jacobian_matrix = open('../inputFiles/jacobian_chaos_cr.txt', 'w')
num_species = ts.shape[1]
for i in range(0,ts.shape[0]):
	if i%200 == 0:
		f_jacob = nd.Jacobian(dX_dt)(np.squeeze(np.asarray(ts[i,:])))
		g.write('%f %f %f %lf %lf\n' % (ts[i, 0], ts[i, 1], ts[i, 2],ts[i,3], ts[i,4]))
		for u in range(0,num_species):
			for z in range(0,num_species):
			    	jacobian_matrix.write('%lf ' % (f_jacob[u,z]))
		jacobian_matrix.write('\n')
g.close()
jacobian_matrix.close()
#### Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.plot(ts[:,0], ts[:,1], ts[:,2], color = 'b')
fig = plt.figure()
plt.plot(first_trace)
plt.show()

