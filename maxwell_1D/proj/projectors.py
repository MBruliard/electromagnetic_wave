__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-12"
__aim__ = "Projectors used into DeRham sequence for maxwell 1D"


# --- Importing modules --- #
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve


from mat1d.mass import *
from space.spaces import coeffDspan

from spl.fem.splines import SplineSpace 
from spl.utilities.quadratures import gauss_legendre
from spl.linalg.stencil        import StencilVector, StencilMatrix



class FirstProjector:
	"""
		First Projector : allows to approach the ordinates corresponding abscissae of greville
		
		Attributs:
		----------
			_values: np.array 
			
			_size: number of values
		
	"""
	def __init__ (self, space,  func, coef=False, t=0.):	
		"""
		
		"""
		
		# ...
		M = MassMatrix(space, coef=coef)
		M = M.tocsc()
		
		size = M.shape[0]
		
		[s1] = space.vector_space.starts
		[e1] = space.vector_space.ends
		[p1] = space.vector_space.pads
		
		k1 = space.quad_order
		spans1 = space.spans
		weights1 = space.quad_weights
		basis1 = space.quad_basis
		points1 = space.quad_points
		
		
		b = np.zeros(size)
		
		for ie1 in range(s1, e1-p1 + 1):
			i_span_1 = spans1[ie1]
			
			for il_1 in range(0, p1+1):
				il1 = i_span_1 - p1 + il_1
				
				for g1 in range(0, k1):
					
					if (coef):
						ni = basis1[ie1, il_1, 0, g1] * coeffDspan(space, il1)
					else:
						ni = basis1[ie1, il_1, 0, g1]
					
					wvol = weights1[ie1, g1]
					pt = points1[ie1, g1]
					
					b[il1] = b[il1] + (wvol*ni*func(pt, t))
					
		
		etild = spsolve(M, b)
		# ...

		self._values = etild
		
	# ...
	
	def tofield (self):
		# TODO
		print("Not implemented yet")
	# ...
	
	def tondarray(self):
		return self._values
	
	def shape (self):
		return self._values.shape
	# ...





