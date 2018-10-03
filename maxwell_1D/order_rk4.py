__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-07-27"


# ----- modules importating
import numpy as np
from math import pi
import matplotlib.pyplot as plt


from mat1d.mass import *
from mat1d.rotational import *
from mat1d.discretdev import *

from scheme.eulerexpl import *
from scheme.rk4 import *
from proj.projectors import *

from userfunc.norm import *



# ----- MAIN FUNCTION
if __name__ == "__main__":

	# data
	ne = 20
	grid =  np.linspace(0., 1., ne)
	degree = 3

	
	#epsilon0 = 8.854e-2
	#mu0 = 12.566370e-7
	epsilon0 = 1
	mu0 = 1
	omega = 1.
	
	leftdirichlet = 0.
	rightdirichlet = 0.
	
	Eanalytic = lambda x, t: np.sin(2*np.pi*omega*x)*np.cos(omega*t)
	Banalytic = lambda x, t: np.cos(2*np.pi*omega*x)*np.sin(omega*t)
	
	# creating spaces
	V, W = creationSpaces(degree, grid)
	
	greville  = V.greville


	# creating startE and startB
	startE = FirstProjector(V, Eanalytic)
	startE = startE.tondarray()
	startB = FirstProjector(W, Banalytic, coef=True)
	startB = startB.tondarray()
	
	exple = startE
	explb = startB

	
	# - Time evolution
	deltaT = [1e-6, 1e-5, 1e-4]
	final_time = 0.05
	err = []
	
	
	# compute exact solution
	efieldex = FirstProjector(V, Eanalytic, t=final_time)
	efieldex = efieldex.tondarray()
	bfieldex = FirstProjector(W, Banalytic, coef=True, t=final_time)
	bfieldex = bfieldex.tondarray()
	
	
	for dt in deltaT:
		t = 0.

		# creating startE and startB
		startE = FirstProjector(V, Eanalytic)
		startE = startE.tondarray()
		startB = FirstProjector(W, Banalytic, coef=True)
		startB = startB.tondarray()
		
			
		# compute euler explicit
		rkfour = SchemeRKfour1D(V, W, startE, startB,  dt, final_time, epsilon0, mu0, dirichlet=True)
		rkfoure = rkfour.efield()
		rkfourb = rkfour.bfield()
		
		
		
		
		# -- compute the error of the scheme
		quad_values_E = quad_values_from_vector(V, rkfoure)
		quad_values_Eexact = quad_values_from_vector(V, efieldex)
		error = normQuadratic(quad_values_E, quad_values_Eexact, V)
		
		print("dt= ", dt)
		print("error L2 ",  error)
		err.append(error)
	# ...
	
	
	# --  creating the error graphic
	order1 = []
	order2 = []
	order3 = []
	for i in range(len(deltaT)):
		order1.append(deltaT[i] )
		order2.append(deltaT[i]**2)
		order3.append(deltaT[i]**3)
	
	plt.loglog()
	plt.plot(deltaT, order1, label='order 1', color='k')
	plt.plot(deltaT, order2, label='order 2', color='b')
#	plt.plot(deltaT, order3, label='order 3', color='c')
	plt.plot(deltaT, err, label='error curve L2', color='r')
	plt.legend()
	plt.title('Error curve of electrical field with RK4\nB-Splines degree: '+str(degree))
	plt.ylabel('log(Error)')
	plt.show()

#	dt = 1e-6
#	# -- creating E vector at time t0
#	oldE = projector_to_l2(V, Eanalytic)
#	oldE = listToVectorColumn(oldE)

#	# -- creating B vector at time t0
#	oldB = projector_to_l2(W, Banalytic)
#	oldB = listToVectorColumn(oldB)
#	
#	t = 0.
#	
#	while (t < final_time):
#			
#		E_eval, B_eval = euler_exp(oldE, oldB, dt, Mdirichlet, Rdirichlet, Gdirichlet, mu0, epsilon0)

#		t = t +dt		
#		oldE = E_eval
#		oldB = B_eval
#		
#		
#		# --  calcul de l'erreur 
#		quad_values_Eexact = quad_values_from_function(V, Eanalytic, t)
#		quad_values_E = quad_values_from_vector(V, oldE)
#		norme = normQuadratic(quad_values_Eexact, quad_values_E, V)
#		
##		print(t, norme)
#		
#			
#	# -- end while
#	print('dt= ', dt)
#	print("norm= ", norme)




	
