__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-12"
__aim__ = "Projectors used into DeRham sequence for maxwell 1D"


# --- Importing modules --- #
import numpy as np
from scipy.sparse import csc_matrix


from mat1d.mass import *
from space.spaces import getNumberOfSplines


from spl.fem.splines import SplineSpace 
from spl.utilities.quadratures import gauss_legendre
from spl.linalg.stencil        import StencilVector, StencilMatrix



class InterpolationMatrix:
	"""
	
	"""

	def __init__ (self, V, coeff=False):
		"""
			We are building the Interpolation Matrix
			
			
		"""
				
			
		# ... sizes
		[s1] = V.vector_space.starts
		[e1] = V.vector_space.ends
		[p1] = V.vector_space.pads
		
		spans_1   = V.spans
		basis_1   = V.quad_basis
		breaks = V.breaks
		
		ne = len(breaks)
		nbS = getNumberOfSplines(V)
		
		# ...
		inter = np.zeros((ne, nbS))
		# ...
		
		for ie1 in range (s1, e1 - p1 + 1):
			print(ie1)
			j_span_1 = spans_1[ie1]
			for jl_1 in range(0, p1+1):
				j = j_span_1 - p1 + jl_1
				
				inter[ie1, j] = basis_1[ie1, jl_1, 0, 0]
		
		for jl_1 in range(0, p1+1):
			inter[-1, j] = basis_1[ie1-1, jl_1, 0, 1]
			
			
		print(inter)
		print(inter.shape)
		
		if (coeff):
			# D_i spans
			# TODO
			print('TODO')
		# ...
		else:
			print("TODO")
			# TODO
		
		
		
		
		# ... 
		self._inter = csc_matrix(inter)
		self._size = inter.shape
	
	# ...
	
	def todense(self):
		
		return self._inter.todense()
	# ...
	
	def dot (self, m):
		return self._inter.dot(m)
	# ...


# ...

class PiGrad:

	def __init__ (V, func, t=0.):
		# TODO
		print('TODO')
	# ...

# ...
