import numpy as np
from wavesolver_2D_ghost import *

####### Plugwave #######

# Ocean floor function
def B(x,y):
	return 0

q_fun = lambda x,y: 0.5

V_fun = lambda x,y: 0

f_fun = lambda x,y,t: 0

def I_fun(x,y):
	if (0.45 < x and x < 0.55):
		val = 1.4
	else:
		val = 1.0

	return val

# Physical constants:
A = 0.5					# amplitude
b = 0					# damping
omega = 1 				# angular frequency 
Lx = 1					# x-length of domain
Ly = 1					# y-length of domain


# Computation parameters
Nx = 100
Ny = 100
dt = 0.001
T = 1.5

dt = (float(Lx)/Nx)/np.sqrt(q_fun(1,1))		#OBS. To do this simulation comment out the stability crit. test in the solver

u = solver_2D(I_fun, V_fun, f_fun, b, q_fun, Lx, Ly, dt, T,Nx ,Ny,B,'none', 'vectorized','animate')



