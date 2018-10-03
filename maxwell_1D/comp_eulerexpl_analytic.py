__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-29"


import numpy as np
from space.spaces import creationSpaces

from util.energy import *
from mat1d.mass import *
from scheme.eulerexpl import *
from scheme.eulerimpl import *
from scheme.rk2 import *
from scheme.rk4 import *
from proj.projectors import *


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
	imple = startE
	implb = startB
	rktwoe = startE
	rktwob = startB
	rkfoure = startE
	rkfourb = startB
	
	deltaT = 5e-5
	maxiter = 1000
	finalT = maxiter* deltaT
	
	t= 0.

	
	for i in range(0, maxiter):
		
		# computing 
		
		# euler explicit
		expl = EulerExplMaxwell1D(V, W, exple, explb,  deltaT, t + deltaT, epsilon0, mu0,  t=t, dirichlet=True)
		exple = expl.efield()
		explb = expl.bfield()
		
		# analytic
		efieldex = FirstProjector(V, Eanalytic, t=t)
		efieldex = efieldex.tondarray()
		bfieldex = FirstProjector(W, Banalytic, coef=True, t=t)
		bfieldex = bfieldex.tondarray()
		
		
		if (i % 200 == 0):
			# ploting results
			plt.plot(greville, exple, label='euler explicit', color='y')
			plt.plot(greville, efieldex, label='analytic', color='r')
			plt.legend()
			plt.title('Electric field\ndt='+str(deltaT)+' t='+str(t))
			plt.savefig('Images/ElectricField1D/EulerExpl/compare_scheme_euler_expl_dt'+str(deltaT)+'_t'+str(t)+'.png')
			plt.clf()
			print("Graph saved in 'Images/ElectricField1D/EulerExpl/compare_scheme_euler_expl_dt"+str(deltaT)+"_t"+str(t)+".png'")
					
		#update
		t = t +deltaT
	# ...
	
