import numpy as np
from wavesolver_2D_ghost import *


###  Verification of constant solution ###

def test_constant(u, u_exact): 
	error = u-u_exact
	error = abs(error)
	error = float(np.sum(error))
	tolerance = 0.0001
	assert error < tolerance
	print 'Test: Constant solution PASSED!'

B = lambda x,y: 0

q_fun = lambda x,y: 1

V_fun = lambda x,y: 0

f_fun = lambda x,y,t: 0

I_fun = lambda x,y: 1.23

#  Parameters:
Ly = 1
Lx = 1
Nx = 100
Ny = 100
dt = 0.001
T = 0.5
b = 0

u = solver_2D(I_fun, V_fun, f_fun, b, q_fun, Lx, Ly, dt, T,Nx ,Ny, B ,'none', 'vectorized','no animate')

u_exact = np.zeros(np.shape(u)) + I_fun(0,0)  # construct exact solution matrix

test_constant(u, u_exact)

