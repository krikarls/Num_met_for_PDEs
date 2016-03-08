
# Exercise 21

from numpy import *
from simulate_pendulum import simulate_pendulum




def eqilimbrium_test():
	x, y, theta, t = simulate_pendulum(0.9, 0, 0, 6, 10, plot=False)
	sum_err = abs(max(x))+abs(max(y))+abs(max(theta))
	tolerance = 1E-5
	if sum_err > tolerance:
		print 'Equilimbrium-test failed'
	else:
		print 'Equilimbrium-test passed'


def vertical_test():
	omega_v = sqrt(beta/(1-beta))
	period = 2*pi/omega_v
	no_periods = 5./omega_v
	timesteps_per_period = omega_v*600
	eps = 0.1

	x, y, theta, t = simulate_pendulum(beta, 0, eps, no_periods, timesteps_per_period, plot=False)

	y_exact = -eps*cos(omega_v*t)

	tolerance = 1E-3
	max_err = max(abs(y-y_exact))

	if max_err > tolerance:
		print 'Vertical-test failed'
	else:
		print 'Vertical-test passed'


beta = 0.99
theta = 9
epsilon = 0
num_periods = 3
time_steps_per_period = 600

eqilimbrium_test()

vertical_test()

simulate_pendulum(beta, theta, epsilon, num_periods, time_steps_per_period, plot=True)