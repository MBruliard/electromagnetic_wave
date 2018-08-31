__author__ = "Margaux BRULIARD"
__date__ = "2018-07-13"
__partners__ = "Ahmed RATNANI", "NMPP"
__purpose__ = "rotational matrix maxwell 1D"

# -- spl packages importing
from spl.fem.splines import SplineSpace 
from spl.utilities.quadratures import gauss_legendre


# -- 
from space.spaces import coeffDspan, getNumberOfSplines

# -- 
import numpy as np
from scipy.sparse import csc_matrix, csr_matrix


class Rotational:

	def __init__ (self, V, W):
		
	
		self._rot = assemblyRotationalMatrix(V, W)
		
		
		self._size = self._rot.shape 
				 
	# ...
	
	
	def dot(self, xs):
		
		return self._rot.dot(xs)
	# ..
	
	def dirichlet(self):
		return self._rot[1:-1, :]
	
	
	
	def todense(self):
		
		return self._rot.todense()
		
	def size(self):
		return self._size


# ...


def assemblyRotationalMatrix(V, W):
	"""
		Assembly the rotationnal matrix between two Spline Spaces @V and @W. 
		
		Entries:
		--------
			V: SplineSpace which represent H^1 space
			
			W: SplineSpace which represent L^2 space 
			
		Returns:
		--------
			R:	sparse matrix
				The rotationnal matrix of size: dim_matrix * (dim_matrix-1)
	"""
	
	# ...
	nbS_v = getNumberOfSplines(V)
	nbS_w = getNumberOfSplines(W)
	R = np.zeros((nbS_v, nbS_w))
	
	
	# ... 
	[s1] = V.vector_space.starts
	[e1] = V.vector_space.ends
	[pv1] = V.vector_space.pads
	[pw1] = W.vector_space.pads

	spans_v = V.spans
	spans_w = W.spans

	basis_v = V.quad_basis
	basis_w = W.quad_basis

	k = V.quad_order
	weight = V.quad_weights

	# --  
	for ie1 in range (0, V.ncells):

		spline_v_index = spans_v[ie1]
		spline_w_index = spans_w[ie1]

		for il in range (0, pv1 + 1):

			i = spline_v_index - pv1 + il

			for jl in range (0, pw1 +1):

				j = spline_w_index - pw1 + jl
				alpha_j = coeffDspan(W, j)
										
				res = 0.0
				for q in range (0, k):
					basis_v_x = basis_v[ie1, il, 1, q]
					basis_w_0 = basis_w[ie1, jl, 0, q]
					
					
					d_w_0 = alpha_j * basis_w_0
					wvol = weight[ie1, q]

					res = res + basis_v_x * d_w_0 * wvol
				R[i, j] = R[i, j] + res
	
	return csr_matrix(R)
# ...


