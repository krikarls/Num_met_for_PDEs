from numpy import *
import matplotlib.pylab as plt

"""
The ODEs to be solved are the following:

x_tt = -(beta/(1-beta))*(1-(beta/L))*x 				# x_tt is the double derivative of x with respect to time
y_tt = -(beta/(1-beta))*(1-(beta/L))*(y-1)-beta 	# these are accelerations

x(0) = (1+epsilon)*sin(theta)
x(0)_t = 0											# initial x-velocity, vx
y(0) = 1 - (1+epsilon)*cos(theta)
y(0)_t = 0											# initial y-velocity, vy 

by setting D = -(beta/(1-beta))*(1-(beta/L)) we obtain a more compact set of equations
L = sqrt(x**2 + (y-1)**2)

x_tt = D*x 				
y_tt = D*(y-1)-beta 	
"""

def simulate_pendulum(beta, theta, epsilon, num_periods, time_steps_per_period, plot=True):

	P = 2*pi
	T = num_periods*P
	dt = P/time_steps_per_period
	omega = 2*pi/P
	theta_0 = theta
	theta = theta*pi/180

	
	t = linspace(0,T, num_periods*time_steps_per_period+1)

	x = zeros(len(t))	# solution arrays
	y = zeros(len(t))

	x[0] = (1+epsilon)*sin(theta)
	y[0] =  1 - (1+epsilon)*cos(theta)

	D = -(beta/(1-beta))*(1-(beta/sqrt(x[0]**2+(y[0]-1)**2)))

	x[1] = 0.5*(dt**2)*D*x[0]+x[0]
	y[1] = 0.5*(dt**2)*D*(y[0]-1)-0.5*beta*dt**2 +y[0]

	for i in range(1,int(len(t)-1)):
		D = -(beta/(1-beta))*(1-(beta/sqrt(x[i]**2+(y[i]-1)**2)))
		x[i+1] = (dt**2)*D*x[i] + 2*x[i] - x[i-1]
		y[i+1] = (dt**2)*D*(y[i]-1) - beta*dt**2 + 2*y[i] - y[i-1]


	theta = arctan(x/(1-y))


	if plot:
		plt.plot(x,y)
		plt.xlabel('x')
		plt.ylabel('y')
		plt.title('Motion of the pendulum')
		plt.axis('equal')
		plt.show()

		plt.plot(t,180*theta/pi)
		plt.xlabel('t')
		plt.ylabel('degrees')
		plt.title('Angular displacement with time')
		plt.show()

		if theta_0 < 10:
			theta_e = theta_0*cos(omega*t)

			plt.plot(t,180*theta/pi,label="Elastic")
			plt.plot(t,theta_e,label="Non-Elastic")
			plt.xlabel('x')
			plt.ylabel('y')
			plt.title('Comparing angular displacement')
			plt.legend()
			plt.show()

	return x, y, theta, t