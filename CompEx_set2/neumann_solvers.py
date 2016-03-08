import numpy as np

# These algorithms store all the u-values stored in a matrix. This is not a ideal, memory wise, but is done due to errors shuffeling arrays. 

# This is a implementation of scheme (50) and (57). HPL. , with Neumann BC
def solver1(I, V, f, q, c_max, L, dt, T, u_e):
	Nt = int(round(T/dt)) 				# Number of last time point(also the number of intervals, since first point has index 0)
	t = np.linspace(0, Nt*dt, Nt+1)   	# Mesh points in time
	dx = (float(c_max)*dt)/0.9 			# Stability criterion for variable wave velocity (0.9 is a safty factor)
	Nx = int(round(L/dx)) 				# Number of last point in space
	x = np.linspace(0, L, Nx+1) 		# Mesh points in space
	dx = x[1]-x[0]

	u = np.zeros([Nt+1,Nx+1])   	# Solution matrix

	K = (dt/dx)**2; DT = dt*dt 		# Constants

	q2 = np.zeros(Nx+1)  			# Computing all q-values(saves FLOPS)
	for i in range(0,Nx+1): 	
		q2[i]=q(x[i])
	q=q2

	import time;  t0 = time.clock()  # Measuring CPU time

	for i in range(0,Nx+1): 	# Set inital condition
		u[0,i] = I(x[i])

	# Computing the the first step using special formulas (only inner spatial points)
	for i in range(1, Nx):
		u[1,i] = u[0,i] + dt*V(x[i]) + 0.25*K*((q[i] + q[i+1])*(u[0,i+1] - u[0,i]) - (q[i] + q[i-1])*(u[0,i] - u[0,i-1])) + 0.5*DT*f(x[i],t[0])

	# Computing the first step boundary values  
	u[1,0] = u[0,0] + dt*V(x[0]) +  0.25*K*((q[0] + q[1])*(u[0,1] - u[0,0]) - (q[0] + q[1])*(u[0,0] - u[0,1])) + 0.5*DT*f(x[0],t[0])
	u[1,Nx] = u[0,Nx] + dt*V(x[Nx]) + 0.25*K*((q[Nx] + q[Nx-1])*(u[0,Nx-1] - u[0,Nx]) - (q[Nx] + q[Nx-1])*(u[0,Nx] - u[0,Nx-1])) + 0.5*DT*f(x[Nx],t[0]) 

	# Computing inner points
	for n in range(2,Nt+1):
		for i in range(1,Nx):
			u[n,i] = -u[n-2,i] + 2*u[n-1,i] + 0.5*K*((q[i] + q[i+1])*(u[n-1,i+1] - u[n-1,i]) - (q[i] + q[i-1])*(u[n-1,i] - u[n-1,i-1])) + DT*f(x[i],t[n-1])

		u[n,0] = -u[n-2,0] + 2*u[n-1,0] + K*0.5*((q[0] + q[1])*(u[n-1,1] - u[n-1,0]) - (q[0] + q[1])*(u[n-1,0] - u[n-1,1])) + DT*f(x[0],t[n-1])
		u[n,Nx] = -u[n-2,Nx] + 2*u[n-1,Nx] + K*0.5*((q[Nx] + q[Nx-1])*(u[n-1,Nx-1]-u[n-1,Nx]) - (q[Nx]+q[Nx-1])*(u[n-1,Nx]-u[n-1,Nx-1])) + DT*f(x[Nx],t[n-1])

	time_cpu = time.clock() - t0 

	u_exact = np.zeros([Nt+1,Nx+1]) 
	if u_e is None or u_e == 0 :
		print 'No exact solution to compare with.'
	else:
		for n in range(0,Nt+1):
			for j in range(0,Nx):
				u_exact[n,j] = u_e(x[j],t[n]) 	# Computing exact solution in all points

	error = abs(u_exact - u)
	mean_err = np.sum(error)/((Nt+1)*(Nx+1)) 	# Mean error per cell

	return u, x, t, time_cpu, mean_err


def solver2(I, V, f, q, c_max, L, dt, T, u_e):
	Nt = int(round(T/dt)) 				# Number of last time point(also the number of intervals, since first point has index 0)
	t = np.linspace(0, Nt*dt, Nt+1)   	# Mesh points in time
	dx = (float(c_max)*dt)/0.9 			# Stability criterion for variable wave velocity (0.9 is a safty factor)
	Nx = int(round(L/dx)) 				# Number of last point in space
	x = np.linspace(0, L, Nx+1) 		# Mesh points in space
	dx = x[1]-x[0]

	u = np.zeros([Nt+1,Nx+1])   	# Solution matrix

	K = (dt/dx)**2; DT = dt*dt 		# Constants

	q2 = np.zeros(Nx+1)  			# Computing all q-values(saves FLOPS)
	for i in range(0,Nx+1): 	
		q2[i]=q(x[i])
	q=q2

	import time;  t0 = time.clock()  # Measuring CPU time

	for i in range(0,Nx+1): 	# Set inital condition
		u[0,i] = I(x[i])

	# Computing the the first step using special formulas (only inner spatial points)
	for i in range(1, Nx):
		u[1,i] = u[0,i] + dt*V(x[i]) + 0.25*K*((q[i] + q[i+1])*(u[0,i+1] - u[0,i]) - (q[i] + q[i-1])*(u[0,i] - u[0,i-1])) + 0.5*DT*f(x[i],t[0])

	# Computing the first step boundary values  
	u[1,0] = u[1,1]
	u[1,Nx] = u[1,Nx-1]

	# Computing inner points
	for n in range(2,Nt+1):
		for i in range(1,Nx):
			u[n,i] = -u[n-2,i] + 2*u[n-1,i] + 0.5*K*((q[i] + q[i+1])*(u[n-1,i+1] - u[n-1,i]) - (q[i] + q[i-1])*(u[n-1,i] - u[n-1,i-1])) + DT*f(x[i],t[n-1])

		u[n,0] = u[n,1]
		u[n,Nx] = u[n,Nx-1]

	time_cpu = time.clock() - t0 

	u_exact = np.zeros([Nt+1,Nx+1]) 
	if u_e is None or u_e == 0 :
		print 'No exact solution to compare with.'
	else:
		for n in range(0,Nt+1):
			for j in range(0,Nx):
				u_exact[n,j] = u_e(x[j],t[n]) 	# Computing exact solution in all points

	error = abs(u_exact - u)
	mean_err = np.sum(error)/((Nt+1)*(Nx+1)) 	# Mean error per cell

	return u, x, t, time_cpu, mean_err

