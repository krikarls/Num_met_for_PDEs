
from dolfin import *
import numpy 
import matplotlib.pyplot as plt

def solver(dt,Nx,T):
	
	mesh = UnitIntervalMesh(Nx)

	V = FunctionSpace(mesh, 'Lagrange', 1)			# Setting up function space

	u0 = Expression('t*(0.5-x[0]/3.0)*x[0]*x[0]', t=0)
	u_1 = interpolate(u0, V) 						# Define initial condition in V

	# define test and trial functions
	u = TrialFunction(V)
	v = TestFunction(V)

	# Variational formulation of the problem:
	# One Picard iteration corresponds to using u from previous time level in alpha
	alpha = lambda u: 1+u**2
	rho = Constant(1.0)
	f = Expression('-rho*pow(x[0],3)/3. + rho*pow(x[0],2)/2. + 8*pow(t,3)*pow(x[0],7)/9. -' + 
                   '28*pow(t,3)*pow(x[0],6)/9. + 7*pow(t,3)*pow(x[0],5)/2. - 5*pow(t,3)*pow(x[0],4)/4.' +
                   ' + 2*t*x[0] - t', rho=rho, t=0)
	a = rho*u*v*dx + inner(dt*alpha(u_1)*nabla_grad(u), nabla_grad(v))*dx
	L = dt*f*v*dx + rho*u_1*v*dx

	u = Function(V)  
	t = dt

	# time step
	while t <= T:
		u0.t=t
		f.t=t

		A = assemble(a)
		b = assemble(L)
		solve(A, u.vector(), b)

		u_1.assign(u)
		t += dt   

	u = u_1
	u_e = Expression('t*(0.5-x[0]/3.0)*x[0]*x[0]', t=t-dt) # to make sure u_e is evaluated at same t
	u_e = interpolate(u_e, V) 	

	e = u.vector().array() - u_e.vector().array()
	E = sqrt(sum(e**2)/u.vector().array().size)

	K = E/dt 
	h = dt

	U = u.vector().array()
	UE = u_e.vector().array()

	return E, h, K, numpy.fliplr([U])[0], numpy.fliplr([UE])[0]


"""      FOR CONVERGENCE RATE
E0=1
h0=1
print 'E', 'h', 'K', 'r'
for n in [4,8,16,32,64]:

	E,h,K, U, U_e = solver((1.0/n)**2,n,1.0)
	print E, h, K , ln(E/E0)/ln(h/h0)
	E0=E
	h0=h
"""

E,h,K, U1, U1_e = solver(0.05,15,1.0)

E,h,K, U3, U3_e = solver(0.05,15,3.0)

E,h,K, U5, U5_e = solver(0.05,15,5.0)

E,h,K, U10, U10_e = solver(0.05,15,10.0)

x = numpy.linspace(0,1,len(U1))

plt.subplot(2,2,1)
plt.plot(x,U1)
plt.plot(x,U1_e)
plt.ylabel('u')
plt.title('T=1')

plt.subplot(2,2,2)
plt.plot(x,U3)
plt.plot(x,U3_e)
plt.ylabel('u')
plt.title('T=3')

plt.subplot(2,2,3)
plt.plot(x,U5)
plt.plot(x,U5_e)
plt.xlabel('x')
plt.ylabel('u')
plt.title('T=5')

plt.subplot(2,2,4)
plt.plot(x,U10)
plt.plot(x,U10_e)
plt.xlabel('x')
plt.ylabel('u')
plt.title('T=10')

plt.xlabel('x')
plt.ylabel('u')
plt.suptitle('Comparing analytical and numerical solution',fontsize=18)
plt.show()




