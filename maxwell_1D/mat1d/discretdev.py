__author__ = "Margaux BRULIARD"
__date__ = "2018-07-13"
__partners__ = "Ahmed RATNANI", "NMPP"
__purpose__ = "Derivative matrix maxwell 1D"

# -- spl packages importing
from spl.fem.splines import SplineSpace 


# -- 
from space.spaces import coeffDspan, getNumberOfSplines

# -- 
import numpy as np
from scipy.sparse import csc_matrix, csr_matrix


class DiscretDev:
	
	def __init__(self, V, W):
		
		self._der = assemblyDerivateMatrix(V)

		
	def dot (self, xs):
		
		return self._der.dot(xs)
	
	def todense(self):
		return self._der.todense()
	
	def size(self):
		return self._size
		
	def dirichlet(self):
		return self._der[:, 1:-1]

# ...


def assemblyDerivateMatrix(V):
	"""
		Entry:
		------
			V: SplineSpace of reference
			
		Output:
		-------
			der: sparse matrix 
	"""	
	
	nbS_v = getNumberOfSplines(V)
	
	M = np.zeros((nbS_v,nbS_v))
	
	for i in range(0, nbS_v):
		M[i,i] = 1.
		# avant c'était 1 ici et -1 plus bas
		if i>0:
			M[i,i-1] = -1.
	# ...
	
	der = M[1:nbS_v, :]
	return csr_matrix(der)

	
# ...
