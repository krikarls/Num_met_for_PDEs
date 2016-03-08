import numpy as np
from wavesolver_2D_ghost import *

####### Standing wave test case #######

# Ocean floor function
def B(x,y):
	return 0

q_fun = lambda x,y: 1

V_fun = lambda x,y: 0

f_fun = lambda x,y,t: 0

u_e = lambda x,y,t : A*np.cos(k_x*x)*np.cos(k_y*y)*np.cos(omega*t)

I_fun = lambda x,y: A*np.cos(k_x*x)*np.cos(k_y*y)

# Physical constants:
A = 1								# amplitude
b = 0								# damping
c_max = 1							# max wave velocity
Lx = 1								# x-length of domain
Ly = 1								# y-length of domain
m_x = 1								# number of x-nodes 
m_y = 1								# number of y-nodes
k_x = (m_x*np.pi)*Lx 
k_y = (m_y*np.pi)/Ly
omega = np.sqrt(k_x**2 + k_y**2)	# angular frequency 

# Computation parameters
Nx = 40
Ny = 40
dt = 0.01
T = 1

p = np.sqrt(2.0/3.0) 	# factors for reducing h = dx*dy*dt
s = 3.0/4.0 				

E = np.zeros(4)
h_vals = np.zeros(4)
h_vals[0] = 0.1 
h_vals[1] = h_vals[0]*0.5
h_vals[2] = h_vals[1]*0.5
h_vals[3] = h_vals[2]*0.5

p = 5

for k in range(1,4):
	Nx = int(round(float(Lx)/(p*h_vals[k-1])))
	Ny = int(round(float(Lx)/(p*h_vals[k-1])))
	dt = h_vals[k-1]
	u, E[k-1] = solver_2D(I_fun, V_fun, f_fun, b, q_fun, Lx, Ly, dt, T,Nx ,Ny,B,u_e, 'vectorized','')
	print 'h= ' ,h_vals[k-1], 'Error= ', E[k-1], 'k=' , k

for i in range(1,3):
	r = np.log(E[i-1]/E[i])/np.log(h_vals[i-1]/h_vals[i])
	print 'Convergence rate: ' , r

