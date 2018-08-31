__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-07-26"
__purpose__ = "Functions for maxwell 1D scheme -- energy of system computing "


import numpy as np

# -- spl packages importing
from spl.fem.splines import SplineSpace 
from spl.utilities.quadratures import gauss_legendre
from spl.linalg.stencil        import StencilVector, StencilMatrix



def computeEnergy (Efield, Bfield, mass, massD , epsilon0, mu0):
	"""
		Compute the energy of the electromagnetic system
		
		Entries:
		--------
			Efield: ndarray - electric field discretized
			
			Bfield: ndarray - magnetic field discretized
			
			mass: sparse matrix -  mass matrix of S^p space
			
			massD: sparse matrix - mass matrix of S^{p-1} space
	"""
	mass = mass.todense()
	massD = massD.todense()
	e = 0.5 * epsilon0 * Efield.dot(mass).dot(Efield)
	b = 0.5 * (1/mu0) * Bfield.dot(massD).dot(Bfield)
	
	Uem = e + b
	
	return  Uem[0, 0]
#...
	




