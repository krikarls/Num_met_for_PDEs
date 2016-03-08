import numpy as np
from wavesolver_2D_ghost import *

####### Manufactured solution  #######


# Ocean floor function
def B(x,y):
	return 0

q_fun = lambda x,y: 4

V_fun = lambda x,y: -np.cos(x)*np.cos(y)

f_fun = lambda x,y,t: 2*(-b*(np.sin(t) + 2*np.cos(t)) + 4*np.sin(t) + 11*np.cos(t))*np.exp(-2*t)*np.cos(x)*np.cos(y)

u_e = lambda x,y,t : (np.cos(t)+np.sin(t))*np.exp(-2*t)*np.cos(x)*np.cos(y)

I_fun = lambda x,y: np.cos(x)*np.cos(y)

# Physical constants:
b = 1					# damping
omega = 1 				# angular frequency 
Lx = 1					# x-length of domain
Ly = 1					# y-length of domain
m_x = 1					# number of x-nodes 
m_y = 1					# number of y-nodes
k_x = (m_x*np.pi)*Lx 
k_y = (m_y*np.pi)/Ly


# Computation parameters
Nx = 100
Ny = 100
dt = 0.001
T = 0.25

u, error = solver_2D(I_fun, V_fun, f_fun, b, q_fun, Lx, Ly, dt, T,Nx ,Ny,B,u_e, 'vectorized',' ')

print error

# well, this is not right.

