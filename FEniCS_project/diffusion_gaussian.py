
from dolfin import *
import numpy as np

# Input parameters
dt = 0.1
T = 5.0
Nx = 20
Ny = Nx
beta = 4.0
sigma = 0.3

#mesh = UnitSquareMesh(Nx, Ny)
mesh =  RectangleMesh(Point(-1,-1), Point(1, 1), 20, 20, "left")

u0 = Expression('exp(-(1/(2*sigma*sigma))*(x[0]*x[0]+x[1]*x[1]))',sigma=sigma)

V = FunctionSpace(mesh, 'Lagrange', 1)		# Setting up function space

u_1 = interpolate(u0, V) 						# Define initial condition in V

# define test and trial functions
u = TrialFunction(V)
v = TestFunction(V)

# Variational formulation of the problem:
# One Picard iteration corresponds to using u from previous time level in alpha
alpha = lambda u: 1+beta*u**2
rho = Constant(50.0)
f = Constant(0)
a = rho*u*v*dx + inner(dt*alpha(u_1)*nabla_grad(u), nabla_grad(v))*dx
L = dt*f*v*dx + rho*u_1*v*dx

u = Function(V)   
t = dt

# time stepping
i=0
while t <= T:
	i += 1
	A = assemble(a)
	b = assemble(L)
	solve(A, u.vector(), b)

	u_1.assign(u)
	t += dt   
u = u_1


#plotter = plot(u)
#plotter.write_png('frame_%04d.png' % i)

file = File("gaussian_diffusion.pvd" )
file << u




