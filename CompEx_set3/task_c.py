from numpy import *
from matplotlib.pyplot import *

c = lambda i: -16.0/(pi*(2*i+1))**3
psi = lambda i,x: sin((2*i+1)*(pi/2.0)*x)

def u(x,N):
	u_accumulate = 0
	for i in range(0,N+1):
		u_accumulate += c(i)*psi(i,x)
	return u_accumulate

x = linspace(0,1,50)
u_exact = 0.5*x**2 - x
u0 = linspace(0,1,50) 	# N = 0
u1 = linspace(0,1,50)	# N = 1
u20 = linspace(0,1,50)	# N = 20

for i in range(0,50):
	u0[i] = u(x[i],0)
	u1[i] = u(x[i],1)
	u20[i] = u(x[i],1)

plot(x,u_exact,label='exact')
plot(x,u0,label='N=0')
plot(x,u1,label='N=1')
plot(x,u20,label='N=20')
title('Approximations using $N$ basis functions $\psi_i$')
legend(loc='upper right')
xlabel('x')
ylabel('u')
show()





