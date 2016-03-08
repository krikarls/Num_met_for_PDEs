
from dolfin import *
import numpy as np

def test_constant(u_array,C):
	error = np.sum(u_array-C) 
	tol = 1e-5
	assert error < tol
	print 'Constant test: PASSED.'

# Input parameters
f = Constant(0);	rho = Constant(1.0); dt = 0.1;T = 1.0;C = 3.0; alpha = lambda u: C
Nx = 10
Ny = 10
Nz = 10

print "Enter dimension:",
dimension = raw_input()

print "Enter 1 for P1 or 2 for P2 elements:",
degree = int(raw_input())

if dimension == '1':
	mesh = UnitIntervalMesh(Nx)
elif dimension == '2':
	mesh = UnitSquareMesh(Nx, Ny)
elif dimension == '3':
	mesh = UnitCubeMesh(Nx, Ny, Nz)
else:
	print 'no valid dimension selected'

u0 = Expression('C', C=C, t=0)

V = FunctionSpace(mesh, 'Lagrange', degree)		# Setting up function space

u_1 = interpolate(u0, V) 						# Define initial condition in V

# define test and trial functions
u = TrialFunction(V)
v = TestFunction(V)

# Variational formulation of the problem:
# One Picard iteration corresponds to using u from previous time level in alpha    
a = rho*u*v*dx + inner(dt*alpha(u_1)*nabla_grad(u), nabla_grad(v))*dx
L = dt*f*v*dx + rho*u_1*v*dx

u = Function(V)   
t = dt

# time stepping
while t <= T:
	u0.t = t

	A = assemble(a)
	b = assemble(L)
	solve(A, u.vector(), b)

	u_1.assign(u)
	t += dt   
u = u_1

u_array = u.vector().array()

test_constant(u_array,C)

plot(u, title='Solution')
interactive()