__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-12"
__aim__ = "Projectors used into DeRham sequence for maxwell 1D"


# --- Importing modules --- #
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve


from mat1d.mass import *
from space.spaces import coeffDspan, getNumberOfSplines

from spl.fem.splines import SplineSpace 
from spl.utilities.quadratures import gauss_legendre
from spl.linalg.stencil        import StencilVector, StencilMatrix



class Histopolation:
	"""
		The histopolation matrix
		
		Attributs:
		----------
			_histo: csc_matrix (sparse matrix)
	"""
	def __init__ (self, V, coef=False, dirichlet=False):
		
		
		# ... sizes
		[s1] = V.vector_space.starts
		[e1] = V.vector_space.ends
		[p1] = V.vector_space.pads
		
		spans_1   = V.spans
		basis_1   = V.quad_basis
		weights_1 = V.quad_weights
		k1        = V.quad_order
		breaks = V.breaks
		
		elem_1 = range(s1, e1 - p1 + 1) 
		ne = len(elem_1)
		nbS = getNumberOfSplines(V)
		
		# ...
		histo = np.zeros((ne, nbS))
		# ...
		
		for ie1 in elem_1:

			j_span_1 = spans_1[ie1]
			
			for jl_1 in range(0, p1+1):
				j = j_span_1 - p1 + jl_1
				
				res = 0.
				for g1 in range(0, k1):
					
					bj_0 = basis_1[ie1, jl_1, 0, g1]
					
					if (coef):
						bj_0 = bj_0 * coeffDspan(V, j)
					
					wol = weights_1[ie1, jl_1]
					
					res = res + (bj_0 * wol)
				# ...
				histo[ie1,j] = histo[ie1, j] + res 
			# ...
		# ...

		if (dirichlet):
			self._histo = csc_matrix(histo[:, 1:-1])
			
		else:		
			self._histo = csc_matrix(histo)
			
	# ...
	
	def solve(self, xs):
		
		return spsolve(self._histo, xs)
	# ...
	
	def shape(self):
		return self._histo.shape
	# ...
	
	def tocsc(self):
		return self._histo
	# ...
	
	def todense(self):
		return self._histo.todense()
	# ...
	
	def dot (self, xs):
		return self._histo.dot(xs)
	# ...
	
# ...


class PiL2:
	"""
		Compute the projector $\pi_{L^2}$
		
		Attributs:
		-----------
			_proj: np.array
			_space: SplineSpace	
	"""
	def __init__ (self, V, func, t=0., coef=False, dirichlet=False, cond=(0., 0.)):
		"""
			
		"""

		if not dirichlet:
			print("PiL2 without dirichlet is not implemented yet")
			exit()
		#..
		
		histo = Histopolation(V, coef=coef, dirichlet=dirichlet)
		
			
		# ... sizes
		[s1] = V.vector_space.starts
		[e1] = V.vector_space.ends
		[p1] = V.vector_space.pads
		
		spans_1   = V.spans
		points_1   = V.quad_points
		weights_1 = V.quad_weights
		k1        = V.quad_order
		breaks = V.breaks
		
		rhs = np.zeros(e1 - p1 +1 - s1)
		
		# TODO: Add the part without Dirichlet Boundary Conditions
		# 		We have to study the one more degree of freedom 
		
		
		# Dirichlet case 
		for ie1 in range(s1, e1 - p1 + 1):

			j_span_1 = spans_1[ie1]
			
			res= 0.
			for g1 in range(0, k1):
			
				pt = points_1[ie1, g1]
				wvol = weights_1[ie1, g1]
				res = res + func(pt, t)*wvol
			# ...
			
			rhs[ie1] = rhs[ie1] + res	
		# ...
		
		
		proj = histo.solve(rhs)
		
		if (dirichlet):
			# we add boundary conditions
			sol = np.zeros(proj.shape[0]+2)
			
			sol[0] = cond[0]
			sol[-1] = cond[-1]
			sol[1:-1] = proj
		
			self._proj = sol
		else:
			self._proj = proj				
		# ...		

		self._space = V

	# ...
	
	def proj (self):
		return self._proj
		
	def shape (self):
		return self._proj.shape
	

# ...
		

