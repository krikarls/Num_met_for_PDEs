# WAVE-PROJECT, solver

#### HOW TO ###
"""
To compute error set u_e as the exact solution as a function(x,y,t), else
set u_e = 'none'.

For animation set visualization = 'animate'.

To visualize with a ocean floor, set B(x,y) to be the desired geometry.
"""
import numpy as np

def solver_2D(I_fun, V_fun, f_fun, b, q_fun, Lx, Ly, dt, T, Nx , Ny, B ,u_e, version, visualization):
	Nt = int(round(T/float(dt))) 		# Number of last time point(also the number of intervals, since first point has index 0)
	t = np.linspace(0, Nt*dt, Nt+1)   	# Mesh points in time
	x = np.linspace(0, Lx, Nx+1) 		# Mesh points in space
	y = np.linspace(0, Ly, Ny+1) 
	dx = x[1]-x[0]						# Making sure dx and dy is the actual distance set by linspace
	dy = y[1]-y[0]

	u = np.zeros([Nx+3,Ny+3])   		# Solution matrix
	u_1 = np.zeros([Nx+3,Ny+3])			# One time-step back
	u_2 = np.zeros([Nx+3,Ny+3])			# Two time-steps back
	u_exact = np.zeros([Nx+1,Ny+1])		# Exact solution 
	ocean_floor = np.zeros([Nx+1,Ny+1])	# Matrix for ocean floor geometry
	Err = np.linspace(0, Nt*dt, Nt+1)	# Array for storing max error of all time-steps					

	q = np.zeros([Nx+3,Ny+3])		# Storing all q-values
	f = np.zeros([Nx+3,Ny+3]) 		# Storing all f-values
	V = np.zeros([Nx+3,Ny+3]) 		# Storing all V-values

	xv = x[:,np.newaxis] 			
	yv = y[np.newaxis,:]

	for j in range(1,len(u[0,:])-1):					# Using for-loop for more flexible 		
		for i in range(1,len(u[:,0])-1):				# input functions I(x,y) and B(x,y).
			u_1[i,j] = I_fun(x[i-1],y[j-1]) 			# Can now contain if-statements.
			ocean_floor[i-1,j-1] = B(x[i-1],y[j-1])

	q[1:-1,1:-1] = q_fun(xv,yv) 	# Pre-compue q values (no time dependence)
	V[1:-1,1:-1] = V_fun(xv,yv)		# Pre-compue V values (no time dependence)
	f[1:-1,1:-1] = f_fun(xv,yv,0)	# Computing f matrix for t=0

	c_max = np.sqrt(np.amax(q))
	test_stability_criterion(c_max,dx,dy,dt) 	# Testing if stability criterion is fulfilled, else: abort

	# Pre-computing constants
	k = (dt*b)/2; Ex = (dt/dx)**2; Ey = (dt/dy)**2; E2x = 0.5*Ex/(1+k); E2y = 0.5*Ey/(1+k); D = (1-k)/(1+k)

	u_1,V,f,q = set_ghost_points(u_1,V,f,q,Nx,Ny,'initiate')	# Set ghost points	

	if visualization == 'animate':
		animator(x,y,t,u_1,0,ocean_floor)

	if version == 'vectorized':		################# Vectorized ################# 
		print 'Vector solver chosen'
		
		# Computing first step: t[1]
		u[1:-1,1:-1] = 0.25*Ex*( (q[1:-1,1:-1]+q[2:,1:-1])*(u_1[2:,1:-1]-u_1[1:-1,1:-1]) - (q[0:-2,1:-1]+q[1:-1,1:-1])*(u_1[1:-1,1:-1]-u_1[0:-2,1:-1]) ) + \
						0.25*Ey*( (q[1:-1,1:-1]+q[1:-1,2:] )*(u_1[1:-1,2:]-u_1[1:-1,1:-1]) - (q[1:-1,0:-2]+q[1:-1,1:-1])*(u_1[1:-1,1:-1]-u_1[1:-1,0:-2]) ) \
						+ u_1[1:-1,1:-1] + (1-k)*dt*V[1:-1,1:-1] + 0.5*(dt**2)*f[1:-1,1:-1] 
		
		f[1:-1,1:-1] = f_fun(xv,yv,t[1]) 			# Update f
		u, f = set_ghost_points(u,V,f,q,Nx,Ny,' ')	# Set ghost points for t[1]
		u_2[:] = u_1; u_1[:] = u					# Shuffle arrays
		
		if visualization == 'animate':
			animator(x,y,t,u,1,ocean_floor) 

		# Computing next steps: t[n+1]
		for n in range(2,Nt+1):
			
			u[1:-1,1:-1] = E2x*( (q[1:-1,1:-1]+q[2:,1:-1])*(u_1[2:,1:-1]-u_1[1:-1,1:-1]) - (q[0:-2,1:-1]+q[1:-1,1:-1])*(u_1[1:-1,1:-1]-u_1[0:-2,1:-1]) ) \
							+ E2y*( (q[1:-1,1:-1]+q[1:-1,2:] )*(u_1[1:-1,2:]-u_1[1:-1,1:-1]) - (q[1:-1,0:-2]+q[1:-1,1:-1])*(u_1[1:-1,1:-1]-u_1[1:-1,0:-2]) )  \
							+ (2/(1+k))*u_1[1:-1,1:-1] - D*u_2[1:-1,1:-1] +((dt**2)/(1+k))*f[1:-1,1:-1]

			f[1:-1,1:-1] = f_fun(xv,yv,t[n])			# Compute f(n) and set ghost points 
			u, f = set_ghost_points(u,V,f,q,Nx,Ny,' ')	# before new time-step 
			u_2[:] = u_1; u_1[:] = u					# Shuffle arrays

			if visualization == 'animate':
				animator(x,y,t,u,n,ocean_floor)
			
			if u_e == 'none':
				a=0 	# a stupid "do nothing argument"
			else:
				u_exact = u_e(xv,yv,t[n])
				error = abs(u_exact-u[1:-1,1:-1])
				Err[n] = np.amax(error) 				# max error per time-step
		
		u = u[1:-1,1:-1]

		if u_e == 'none': 		
			return u
		else:
			return u, np.amax(Err) 

	
	elif version == 'scalar': 		################### Scalar ################### 
		print 'Scalar solver chosen' 
		
		# Computing first step: t[1]
		for j in range(1,len(u[0,:])-1):
			for i in range(1,len(u[:,0])-1):

				u[i,j] = E2x*( (q[i,j]+q[i+1,j])*(u_1[i+1,j]-u_1[i,j])-(q[i-1,j]+q[i,j])*(u_1[i,j]-u_1[i-1,j]) ) + \
							E2y*( (q[i,j]+q[i,j+1])*(u_1[i,j+1]-u_1[i,j]) - (q[i,j-1]+q[i,j])*(u_1[i,j-1]-u_1[i,j]) ) + \
							u_1[i,j] + (1-k)*dt*V[i,j] + ((dt**2)/2.0)*f[i,j]

		f[1:-1,1:-1] = f_fun(xv,yv,t[1])			# Compute f(n) and set ghost points 
		u, f = set_ghost_points(u,V,f,q,Nx,Ny,' ')	# before new time-step 
		u_2[:] = u_1; u_1[:] = u					# Shuffle arrays

		if visualization == 'animate':
			animator(x,y,t,u,1,ocean_floor) 


		# Computing next steps: t[n+1]
		for n in range(2,Nt+1):
			for j in range(1,len(u[0,:])-1):
				for i in range(1,len(u[:,0])-1):

					u[i,j] = E2x*( (q[i,j]+q[i+1,j])*(u_1[i+1,j]-u_1[i,j])-(q[i-1,j]+q[i,j])*(u_1[i,j]-u_1[i-1,j]) ) + \
							E2y*( (q[i,j]+q[i,j+1])*(u_1[i,j+1]-u_1[i,j]) - (q[i,j-1]+q[i,j])*(u_1[i,j-1]-u_1[i,j]) ) + \
							(2.0/(1+k))*u_1[i,j] - D*u_2[i,j] + ((dt**2)/(1+k))*f[i,j]

			f[1:-1,1:-1] = f_fun(xv,yv,t[n])			# Compute f(n) and set ghost points 
			u, f = set_ghost_points(u,V,f,q,Nx,Ny,' ')	# before new time-step 
			u_2[:] = u_1; u_1[:] = u					# Shuffle arrays

			if visualization == 'animate':
				animator(x,y,t,u,n,ocean_floor)


		if test_case == 'constant':
			test_constant(u[1:-1,1:-1],I_fun)

		return u[1:-1,1:-1] 

	else:							
		print 'No valid solver option chosen'

######################## END SOLVER FUNCTION ########################



def set_ghost_points(u,V,f,q,Nx,Ny,step): 		

	u[0,1:-1] = u[2,1:-1]
	u[Nx+2,1:-1] = u[Nx,1:-1]
	u[1:-1,0] = u[1:-1,2]
	u[1:-1,Ny+2] = u[1:-1,Ny] 

	f[0,1:-1] = f[2,1:-1]
	f[Nx+2,1:-1] = f[Nx,1:-1]
	f[1:-1,0] = f[1:-1,2]
	f[1:-1,Ny+2] = f[1:-1,Ny] 

	if step == 'initiate':			# Since V and q has no time dependence they only
		V[0,1:-1] = V[2,1:-1]		# to be set in the begining
		V[Nx+2,1:-1] = V[Nx,1:-1]
		V[1:-1,0] = V[1:-1,2]
		V[1:-1,Ny+2] = V[1:-1,Ny] 

		q[0,1:-1] = q[2,1:-1]
		q[Nx+2,1:-1] = q[Nx,1:-1]
		q[1:-1,0] = q[1:-1,2]
		q[1:-1,Ny+2] = q[1:-1,Ny] 

		return u,V,f,q

	else:
		return u, f

def test_stability_criterion(c_max,dx,dy,dt):
	stability_factor = (1/c_max)*(1.0/np.sqrt((1/dx**2)+(1/dy**2)))
	assert  dt < stability_factor
	print 'Stability criterion fulfilled.'

def animator(x,y,t,u,n,ocean_floor):
	if n % 10 == 0:

		from mpl_toolkits.mplot3d import Axes3D
		from matplotlib import cm
		from matplotlib.ticker import LinearLocator, FormatStrFormatter
		import matplotlib.pyplot as plt

		fig = plt.figure()
		ax = fig.gca(projection='3d')
		X, Y = np.meshgrid(y, x)
		Z = u[1:-1,1:-1]				# Water surface
		Z2 = ocean_floor				# Ocean floor

		# Different plot styles
		ax.set_zlim(0, 1.5)
		surf = ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
		surf = ax.plot_surface(X, Y, Z2, color='0.75')
		# surf = ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
		# surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.ocean, linewidth=0, antialiased=False)

		plt.title(['t=%f' % t[n]])
		plt.savefig('/home/krister/Documents/INF5620/wave_project/movie_frames/frame_%04d.png' % n)

		# TO MAKE MOVIE: 
		# run: "convert -delay 10 -loop 0 frame_*.png animation.gif" from terminal 
