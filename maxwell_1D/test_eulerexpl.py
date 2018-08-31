__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-29"


import numpy as np
from space.spaces import creationSpaces

from scheme.eulerexpl import *
from proj.projectors import *


if __name__ == "__main__":

	# data
	ne = 25
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
	
	# creating startE and startB
	startE = FirstProjector(V, Eanalytic)
	startE = startE.tondarray()
	startB = FirstProjector(W, Banalytic, coef=True)
	startB = startB.tondarray()
	
	
	deltaT = 5e-5
	finalT = 0.1
	
	# computing 
	res = EulerExplMaxwell1D(V, W, startE, startB,  deltaT, finalT, epsilon0, mu0,  t=0., savefig=True, dirichlet=True)
	
	
	mass = MassMatrix
#	print("Field E: ")
#	print(res.efield())
#	print("Field B: ")
#	print(res.bfield())
#	
	efield = res.efield()
	bfield = res.bfield()
	
	efieldex = FirstProjector(V, Eanalytic, t=finalT)
	efieldex = efieldex.tondarray()
	bfieldex = FirstProjector(W, Banalytic, coef=True, t=finalT)
	bfieldex = bfieldex.tondarray()

	
	
	# ploting ending 
	plt.plot(efield, label='approximate solution')
	plt.plot(efieldex, label='exact solution')
	plt.plot(startE, label='initial field')
	plt.legend()
	plt.title('Euler Explicit Scheme Electric field\ndt='+str(deltaT)+' t='+str(finalT))
	plt.show()


	
# ...
