__author__ = "Margaux BRULIARD"
__date__ = "22.08.2018"
__partners__ = "Ahmed RATNANI", "NMPP", "IPP"


import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import inv
from matplotlib import pyplot as plt

from space.spaces import *
from mat1d.mass import MassMatrix
from mat1d.rotational import Rotational
from mat1d.discretdev import DiscretDev

from util.vector import *

class SchemeRKfour1D :
	
	def __init__ (self, V, W, startE, startB,  deltaT, finalT, epsilon0, mu0,  t=0., savefig=False, dirichlet=False, cond=(0., 0.)):
		"""
			Entries:
			--------
				V, W : SplineSpace 
				
				startE: ndarray (n, 1) - first field E
				
				startB: ndarray(m, 1) - first field B
				
				finalT: float - final time
				
				deltaT:  time step
				
				savefig: boolean to decide if we save the differents graphics
				
				dirichlet: boolean 
				
				cond: tuple - cond[0] = left dirichlet condition, cond[1] = right dirichlet condition 
			
			Output:
			-------
				
		"""
		mass = MassMatrix(V)
		invmass = mass.inv()
		rot = Rotational(V, W)
		dis = DiscretDev(V, W)
		
		if (dirichlet==True):
			
			mass = mass.dirichlet()
			rot = rot.dirichlet()
			dis = dis.dirichlet()
			
			invmass = inv(mass)
			
			startE = startE[1:-1]
			
		# ...
		
		nbiter = int(finalT/deltaT)
		i = 0
		
		startE = listToVectorColumn(startE)
		startB = listToVectorColumn(startB)
		
		while (t < finalT):
		
			# TODO: a revoir le signe -> selon la construction de rot 
			ek1 = deltaT * derE(invmass, rot, startB, mu0, epsilon0)
			bk1 = deltaT* derB(dis, startE)
			
			ek2 = deltaT * derE(invmass, rot, startB + 0.5*bk1, mu0, epsilon0)
			bk2 = deltaT * derB(dis, startE + 1/2 * ek1)
			
			ek3 = deltaT * derE(invmass, rot, startB + 1/2 * bk2, mu0, epsilon0)
			bk3 = deltaT* derB(dis, startE + 1/2*ek2)
			
			ek4 = deltaT* derE(invmass, rot, startB +  bk3, mu0, epsilon0)
			bk4 = deltaT * derB(dis, startE + ek3)
			
			
			
#			E = startE + deltaT/6 * (ek1 + 2*ek2 + 2* ek3  + ek4)
#			B = startB + deltaT/6 * (bk1 + 2*bk2 + 2*bk3 + bk4)
			E = startE + 1/6 * (ek1 + 2*ek2 + 2* ek3  + ek4)
			B = startB + 1/6 * (bk1 + 2*bk2 + 2*bk3 + bk4)
			
			# savefig == True -> we save in some graphs the evolvment of E and B
			if (savefig):
				if (i % 100 == 0):
					
					appearE = E
					if (dirichlet):
						appearE = addDirichletConditions(appearE, cond[0], cond[-1])
					
					plt.plot(V.greville, appearE, label='field E', color='g')
					plt.plot(W.greville, B, label='field B', color='y')
					plt.legend()
					plt.xlabel('$\Omega$')
					plt.ylim(-1.3, 1.3)
					plt.title('Approximation of the electromagnetic fields with RK4 scheme\nt='+str(t)+' dt='+str(deltaT))
					plt.savefig('Images/rk4/rk4_dt'+str(deltaT)+'_it'+str(int(i/100)+1)+'.png')
					print("Image saved in 'Images/rk4/rk4_dt"+str(deltaT)+"_it"+str(int(i/100)+1)+".png'") 
					plt.clf()
					
			# update
			startE = E
			startB = B
			t = t + deltaT
			i= i +1
			
		# ...

		if (dirichlet):
			E = addDirichletConditions(E, cond[0], cond[-1])
		
		
		self._efield = E
		self._bfield = B
	# ...

	def efield(self):
		return self._efield
	
	def bfield(self):
		return self._bfield
	
# ...


def derE (invmass, rot, bfield, mu0, epsilon0):
	
	return 1/(mu0 * epsilon0) * invmass.dot(rot).dot(bfield)

# ...

def derB (dis, efield):

	return dis * efield
#...


	
