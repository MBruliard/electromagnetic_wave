__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-22"
__purpose__ = "Kronecker product between two matrices"


import numpy as np
from util.mattransform import fourToTwoD, nonameyet


class Kronecker:

	def __init__ (self, A, B, twodim=False):
		"""
			Compute Kronecker Product between A et B
			
			Requires:
			---------
				@A ndarray 2D
				@B ndarray 2D
				@2D boolean: 'false'= ndarray 4D, 'true' = ndarray 2D
			
			Returns:
			--------
				@K kronecker product ndarray 2D or 4D
		"""

		
		(m1, m2) = A.shape
		(n1, n2) = B.shape
		K = np.zeros((m1, m2, n1, n2))
		
		for i in range(0, m1):
			for j in range(0, m2):
				K[i, j, :, :] = A[i, j]*B
		
		if(twodim):
			K = fourToTwoD(K)	
		
		self.Kron = K
	
	def matrix(self):
		return self.Kron
