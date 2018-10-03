__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-20"
__purpose__ = "discrete derivation matrix (G) in Maxwell 2D"

import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, kron

from space.spaces import *
from util.kronecker import *
from fieldtest import *

from spl.fem.splines import SplineSpace
from spl.linalg.stencil import StencilVector, StencilMatrix
from spl.fem.tensor  import TensorFemSpace


class DiscreteDev:


	def __init__ (self, V):
		
		[V1, V2] = [W for W in V.spaces]
		v1 = getNumberSplineInSpace(V1)
		v2 = getNumberSplineInSpace(V2)
		

		G1 = kernelG1(v1, v2)
		G2 = kernelG2(v2, v1)		
		
#		G1 = Kronecker(Dis, Id, twodim='false')
#		G2 = Kronecker(-Id, Dis, twodim='true')
#		
		
		self._dis = [G1, G2]
		self.size = (G1.shape, G2.shape)
		
	def dot (self, x):
		"""
			Compute G*x
			
			Entry:
			------
				x: Field 1D
				
			Output:
			-------
				y: Field 2D
		"""
			
		y1 = self._dis[0].dot(x.dir())
		y2 = self._dis[1].dot(x.dir())

#		y1 = self._dis[0].dot(x._values)
#		y2 = self._dis[-1].dot(x._values)
#		
#		y1 = nonameyet(self._dis[0], x._values)
#		y2 = nonameyet(self._dis[-1], x._values)
		
		return Field2d(y1, y2)
	# ...
	
#	
#	def todense(self):
#		
#		return [self.G[0].todense(), self.G[1].todense()]
	
# ...




def assembly_discrete_derivation_1d (n):
	"""
		Compute the discrete derivation matrix 1D case
		
		Entry:
		------
			n: INTEGER - number of column of G matrix
			
		Output:
		-------
			G: sparse matrix 2D of size (n - 1, n)
	"""
	
	G = np.zeros((n,n))
	for i in range(0, n):
		G[i,i] = 1.
		
		if i>0:
			G[i,i-1] = -1.
	
	return G[1:n,:]
	
	
	
def kernelG1 (sizeDis, sizeId):
	"""
		Compute the matrix G1 
		
		Entry:
		------
			sizeId: INTEGER 
			
		
		Output:
		-------
			G1: sparse matrix $ G1 = Dis \otimes Id $ 
	"""
	Dis = assembly_discrete_derivation_1d(sizeDis)
	Id =  np.eye(sizeId)
	
	G1 = kron(Dis, Id)
#	G1 = Kronecker(Dis, Id).matrix()
	
	return G1
	

def kernelG2 (sizeDis, sizeId):
	"""
		Compute the matrix G1 
		
		Entry:
		------
			sizeId: INTEGER 
			
		
		Output:
		-------
			G2: sparse matrix $ G1 = -Id \otimes Dis $ 
	"""
	
	Dis= assembly_discrete_derivation_1d(sizeDis)
	Id = np.eye(sizeId)
	
	G2 = kron(-Id, Dis)
#	G2 = Kronecker(-Id, Dis).matrix()
	
	return G2
	
	
	
	
	
	
	
