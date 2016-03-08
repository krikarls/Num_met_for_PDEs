import numpy as np
import scitools.std as plt
import sympy as sym
from neumann_solvers import *

PI, X, t_, l = sym.symbols('PI, X, t_, l')

# Function for fitting source t
def source_term(Q,U):
	F = sym.diff(U,t_,t_) - sym.diff(Q,X)*sym.diff(U,X) - Q*sym.diff(U,X,X)
	return sym.simplify(F)

####### Exercise b ####### 

U = sym.cos((PI*X) / l)*sym.cos(t_)
Q = 1+sym.cos(PI*X/l)

F = source_term(Q,U) # Symbolic source term

def f(x,t): 	# Function source term
	return (-2*np.pi**2*np.sin(np.pi*x/L)**2 + np.pi**2*np.cos(np.pi*x/L) + np.pi**2 - L**2*np.cos(np.pi*x/L))*np.cos(t)/L**2

V = lambda x: 0
I = lambda x: np.cos(np.pi*x/L)
q = lambda x: 1+np.cos(np.pi*x/L)				
u_e = lambda x,t: np.cos(np.pi*x/L)*np.cos(t)

L = 1
c_max = np.sqrt(2)
dt = 0.2
T = 2

for k in range(0,5):
	u, x, t, cpu_time, err = solver1(I, V, f, q, c_max, L, dt*(0.5)**k, T, u_e)
	print 'dt = ', dt*(0.5)**k , '   Mean cell error = ', err


# For movie making
"""
u, x, t, cpu_time, err = solver1(I, V, f, q, c_max, L, dt, T, u_e)
u_ex = np.zeros([len(t),len(x)])

for n in range(0,len(t)):
	for i in range(0,len(x)):
		u_ex[n,i] = u_e(x[i],t[n])

for h in range(0, len(u[:,1])):
	plt.plot(x, u[h,:],x,u_ex[h,:],legend=("Numerical","Exact"), xlabel='x', ylabel='u',axis=[0, L, -2, 2],title='t=%f' % t[h])
	plt.savefig('frame_%04d.png' % int(h))

"""