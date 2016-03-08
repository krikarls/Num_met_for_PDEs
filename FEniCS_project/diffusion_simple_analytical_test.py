
from dolfin import *
import numpy 


def solver(dt,Nx):
	# parameters
	Ny = Nx; T = 1.0
	
	mesh = UnitSquareMesh(Nx, Ny)

	V = FunctionSpace(mesh, 'Lagrange', 1)			# Setting up function space

	u0 = Expression('exp(-pi*pi*t)*cos(pi*x[0])', t=0)
	u_1 = interpolate(u0, V) 						# Define initial condition in V

	# define test and trial functions
	u = TrialFunction(V)
	v = TestFunction(V)

	# Variational formulation of the problem:
	# One Picard iteration corresponds to using u from previous time level in alpha
	alpha = lambda u: 1
	f = Constant(0) 
	rho = Constant(1.0)
	a = rho*u*v*dx + inner(dt*alpha(u_1)*nabla_grad(u), nabla_grad(v))*dx
	L = dt*f*v*dx + rho*u_1*v*dx

	u = Function(V)  
	t = dt

	# time step
	while t <= T:
		u0.t=t

		A = assemble(a)
		b = assemble(L)
		solve(A, u.vector(), b)

		u_1.assign(u)
		t += dt   

	u = u_1
	u_e = Expression('exp(-pi*pi*t)*cos(pi*x[0])', t=t-dt) # to make sure u_e is evaluated at same t
	u_e = interpolate(u_e, V) 	

	e = u.vector().array() - u_e.vector().array()
	E = sqrt(sum(e**2)/u.vector().array().size)

	K = E/dt 
	h = dt

	return E, h, K

E0=1.0
h0=1.0
print 'E', 'h', 'K', 'r'
for n in [4,8,16,32,64]:
	E,h,K = solver((1.0/n)**2,n)
	if n != 4: 		# need two E/h pairs before convergence rate makes sense
		print E, h, K , ln(E/E0)/ln(h/h0)
	E0=E
	h0=h

"""
E,h,K = solver(0.05,10)

print E
print h
print K

E,h,K = solver(0.025,10)

print E
print h
print K

print 'E', 'h', 'K'
E,h,K = solver(0.05,4)
print E, h, K
E,h,K = solver(0.025,5)
print E, h, K
E,h,K = solver(0.0125,9)
print E, h, K
E,h,K = solver(0.00625,13)
print E, h, K
E,h,K = solver(0.003125,18)
print E, h, K
"""

"""
h N
0.05 4
0.025 6
0.0125 9
0.00625 13
0.003125 18
"""


