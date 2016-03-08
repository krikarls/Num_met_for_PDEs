

import matplotlib.pyplot as plt
from numpy import *
import sympy as sym

def vib_solver(dt, I, V, w, T, b, plot=False):

	# Source-term
	# Linear if b=0, else quadratic
	def f(n):
		return (2*b**2) + (w**2)*(b*(n*dt)**2+V*(n*dt)+I)

	N = float(T/dt)+1 	# number of mesh-points
	t = linspace(0,T,N)
	u = zeros(N)

	u[0] = I
	u[1] = u[0] +V*dt +0.5*(dt**2)*(f(0)-(w**2)*u[0])

	for i in range(1,int(N-1)):
		u[i+1] = f(i)*dt**2 - (dt**2)*(w**2)*u[i] + 2*u[i] - u[i-1]

	u_e = b*t**2 +V*t + I	

	if plot:
		plt.plot(t,u_e,label="Exact")
		plt.plot(t,u,label="Numerial solution")
		plt.xlabel('t')
		plt.ylabel('u')
		plt.title('Comparing numerical and exact solution for dt=%.3f' % dt)
		plt.legend()
		plt.show()

	if plot == False:
		return u, u_e, t