import numpy as np
from wavesolver_2D_ghost import *

####### Wave propagation over uneven ocean floor  #######

# function defning the ocean floor B(x,y)

def B(x,y):			# Dramatic Guassian hill		
	return 0.7*np.exp(- ((x-0.5)/0.1)**2 -  ((y - 0.5)/0.1)**2 )

def q_fun(x,y):
	return g*(1-B(x,y))

V_fun = lambda x,y: 0

f_fun = lambda x,y,t: 0

u_e = lambda x,y,t : A*np.cos(k_x*x)*np.cos(k_y*y)*np.cos(omega*t)

I_fun = lambda x,y: 1 + np.exp(-(x/0.1)**2)

# Physical constants:
g = 9.81				# gravity
A = 1					# amplitude
b = 0					# damping
c_max = 1				# max wave velocity
omega = 1 				# angular frequency 
Lx = 1					# x-length of domain
Ly = 1					# y-length of domain
m_x = 1					# number of x-nodes 
m_y = 1					# number of y-nodes
k_x = (m_x*np.pi)*Lx 	# wave number
k_y = (m_y*np.pi)/Ly 	# wave number

# Computation parameters
Nx = 100
Ny = 100
dt = 0.001
T = 1

u = solver_2D(I_fun, V_fun, f_fun, b, q_fun, Lx, Ly, dt, T,Nx ,Ny, B ,'none', 'vectorized','animate')


