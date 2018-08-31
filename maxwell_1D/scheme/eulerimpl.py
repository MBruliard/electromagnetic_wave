__author__ = "Margaux BRULIARD"
__date__ = "22.08.2018"
__partners__ = "Ahmed RATNANI", "NMPP", "IPP"


import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import inv, spsolve
from matplotlib import pyplot as plt

import sys

from space.spaces import *
from mat1d.mass import MassMatrix
from mat1d.rotational import Rotational
from mat1d.discretdev import DiscretDev

from util.vector import *


class EulerImplMaxwell1D :
	
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
		
		# .. 
		# TODO:  a revoir la construction du rot et modifier en fonction 
		implmat = assemblyImplicitMatrix(invmass, -rot, dis, deltaT, epsilon0, mu0 )
		
		# TODO check inversibility
#		if (not checkInversibilityImplicitMatrix(invmass, dis, rot, deltaT, epsilon0, mu0)):
#			print("ERROR: implicit matrix is singular")
#			sys.exit()
			
		# .. build vector of computing:
		n = startE.shape[0]
		b = group_two_vectors(startE, startB)
		# ..	
		
		
		i = 0
		
		
		while (t < finalT):
		
			new = spsolve(implmat, b)
			
			# savefig == True -> we save in some graphs the evolvment of E and B
			if (savefig):
				if (i % 100 == 0):
					
					appearE, appearB = ungroup_vector(new, n)
					
					if (dirichlet):
						appearE = addDirichletConditions(appearE, cond[0], cond[-1])
						
							
					plt.plot(appearE, label='field E', color='g')
					plt.plot(appearB, label='field B', color='y')
					plt.legend()
					plt.xlabel('$\Omega$')
					plt.ylim(-1.3, 1.3)
					plt.title('Approximation of the electromagnetic fields with Euler Implicit scheme\nt='+str(t)+' dt='+str(deltaT))
					plt.savefig('Images/EulerImpl/eulerimpl_dt'+str(deltaT)+'_it'+str(int(i/100)+1)+'.png')
					print("Image saved in 'Images/EulerImpl/eulerimpl_dt"+str(deltaT)+"_it"+str(int(i/100)+1)+".png'") 
					plt.clf()
					
			# update
			b = new
			t = t + deltaT
			i= i +1
			
		# ...
		E, B = ungroup_vector(b, n)
		
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


def assemblyImplicitMatrix(invmass, rot, dis, deltaT, epsilon0, mu0 ):
	"""
		Requires:
		---------
			invmass: sparse matrix - inverse of the mass matrix 
			
			rot: Rotational matrix
			
			dis: DiscretDev matrix
			
			deltaT: is a real
			
			epsilon0 and mu0 are constant
		
		Returns:
		--------
			res np.ndarray 2D of size (n+m, n+m)
	"""
		
	# -- define the size of res matrix
	n = invmass.shape[0]	
	m = dis.shape[0]
		
	res = np.zeros((n+m, n+m))
	
	# -- first corner
	first = np.eye(n)
	second = - deltaT/(mu0*epsilon0)*np.dot( invmass.todense(), rot.todense()  )

	third = - deltaT*dis.todense()
	fourth = np.eye(m)


	# -- res matrix
	res[0:n, 0:n] = first
	res[0:n, n:n+m] = second
	
	res[n:n+m, 0:n] = third
	res[n:n+m, n:n+m] = fourth
	 
	 	
	
	return csc_matrix(res);
# ..

def checkInversibilityImplicitMatrix(invM, G, R, deltaT, epsilon0, mu0):
	"""
	
	We know that the implicit matrix is build:
	
		Implicit Matrix = [[ 		Id_n				\frac{\Delta T}[\epsilon_0 \mu_0] M^{-1} R]
								
							- \frac{Delta T}G							Id_{m}					]]
		Require:
		--------
			M: mass matrix in presented in implicit matrix of size (n, n)
			
			R: rotational matrix of size (n, m)
			
			G: discrete derivate matrix of size (m, n)
		
		Return:
		-------
		
			choice: bool 
	"""
	invM = invM.todense()
	
	(n , m) = R.shape
	
	A = np.eye(n)
	D = np.eye(m)
	
	# -- we have to show : D + G A inv(M) R is non-sigular
	alpha = deltaT**2/(epsilon0*mu0)
	Check = D + alpha *  G.dot(invM).dot(R)
	
	
	det = np.linalg.det(Check)
	
	if (det == 0):
		return False
	
	return True
	



