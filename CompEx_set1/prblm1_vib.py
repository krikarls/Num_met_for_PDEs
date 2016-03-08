
# INF5620 - PROBLEM 1 (insert date)
# Krister S. Karlsen

from vib_solver import vib_solver
import sympy as sym
import numpy as np

a,b,V,I,t,dt,w = sym.symbols('a b V I t dt w')

# Un-comment the one you want to use 

# u = V*t + I 						; print 'Now using quadratic'				
u = b*t**2 + V*t + I 				; print 'Now using quadratic'			
# u = a*t**3 + b*t**2 + V*t + I 	; print 'Now using cubic'	

# returns the left hand of the ODE, which is equal to f(t)
def ODE_source_term(u,w):
	return sym.diff(u,t,t)+(w**2)*u

# computes finite difference approx for DtDt
def DtDt_fun(u, dt):	
	return (u.subs(t,t+dt)-2*u+u.subs(t,t-dt))/dt**2	

# computes the residual
def residual_discrete_eq(u,DtDt,f):	
	return DtDt + (w**2)*u - f

# computes the residual from the first step
def residual_first_step(u,V,dt,f,w):
	return u.subs(t,0)+V*dt+(dt**2/2)*(f.subs(t, 0)-(w**2)*u.subs(t,0))- u.subs(t,dt)

DtDt = sym.simplify( DtDt_fun(u, dt) )

f = ODE_source_term(u,w)
print 'The source term is: ', f

res = sym.simplify( residual_discrete_eq(u,DtDt,f) )
print 'The residual is: ', res

first_res = sym.simplify( residual_first_step(u,V,dt,f,w) )
print 'The residual from the first step is: ', first_res


dt = 0.2	; I= 1	; V=1	; w= 6	; T= 5; b=1

vib_solver(dt, I, V, w, T, b, plot=True)


###### Testing convergence ###### 

dt_values = []
E_values = []
for i in range(5):
	u, u_e, t = vib_solver(dt, I, V, w, T, b, plot=False)
	E = np.sqrt(dt*np.sum((u_e-u)**2))
	dt_values.append(dt)
	E_values.append(E)
	dt = dt/2
	r = np.log(E_values[i-1]/E_values[i])/np.log(dt_values[i-1]/dt_values[i])
	print r

###### -- ###### 


