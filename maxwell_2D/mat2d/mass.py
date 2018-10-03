__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-21"
__purpose__ = "assembly mass matrix for Maxwell 2D equation"

from spl.fem.splines import SplineSpace
from spl.linalg.stencil import StencilVector, StencilMatrix
from spl.fem.tensor  import TensorFemSpace

from space.spaces import coefficients_linear_combinaison
from util.mattransform import fourToTwoD

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import inv


class Mass :

	def __init__ (self, V, coeff=False):
	
		if (coeff):
			self._mass = assemblyMassMatrixDspan2D(V)
		else:
			self._mass = assemblyMassMatrix2D(V)
		self._size = self._mass.tocoo().shape
	# ...
	
	
	def tocoo(self):
		return self._mass.tocoo()
	# ...
		
	def tocsc (self):
		return csc_matrix(self._mass.tocoo())
	# ...
	
	def todense(self):
		return self._mass.tocoo().todense()
	# ...
	
	def inv(self):
		
		m = self.tocsc()
		
		return inv(m)
	# ...
	
	
	def dot(self, xs):
		return self._mass.tocsc().dot(xs)
	# ...
	
	def to2Darray(self):
		return csc_matrix(fourToTwoD(self._mass.todense()))
		
	def size(self):
		return self.tocsc().shape

# ...


def assemblyMassMatrix2D(V):
	"""
	 	Require:
	 	--------
	 		@V is a FemTensor V= (V1, V2) with @V1, @V2 are FemSpace
	
		Returns:
		--------
			@M is a StencilMatrix - the mass matrix for V space
	"""

	# ... sizes
	[s1, s2] = V.vector_space.starts
	[e1, e2] = V.vector_space.ends
	[p1, p2] = V.vector_space.pads
	# ...

	# ... seetings
	[k1, k2] = [W.quad_order for W in V.spaces]
	[spans_1, spans_2] = [W.spans for W in V.spaces]
	[basis_1, basis_2] = [W.quad_basis for W in V.spaces]
	[weights_1, weights_2] = [W.quad_weights for W in V.spaces]
	[points_1, points_2] = [W.quad_points for W in V.spaces]

	# ... data structure
	M = StencilMatrix(V.vector_space, V.vector_space)
	# ...

	# ... build matrices
	for ie1 in range(s1, e1+1-p1):
		for ie2 in range(s2, e2+1-p2):
			i_span_1 = spans_1[ie1]
			i_span_2 = spans_2[ie2]
			for il_1 in range(0, p1+1):
				for jl_1 in range(0, p1+1):
					for il_2 in range(0, p2+1):
						for jl_2 in range(0, p2+1):

							i1 = i_span_1 - p1  - 1 + il_1
							j1 = i_span_1 - p1  - 1 + jl_1

							i2 = i_span_2 - p2  - 1 + il_2
							j2 = i_span_2 - p2  - 1 + jl_2

							v_m = 0.0
							for g1 in range(0, k1):
								for g2 in range(0, k2):
									bi_0 = basis_1[ie1, il_1, 0, g1] * basis_2[ie1, il_2, 0, g2]
									bj_0 = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 0, g2]
									
									wvol = weights_1[ie1, g1] * weights_2[ie2, g2]

									v_m = v_m + bi_0 * bj_0 * wvol
							M[i1, i2, j1 - i1, j2 - i2]  += v_m
	# ...
	return M
# ...


# ...
def assemblyMassMatrixDspan2D(V):
	"""
		We compute the mass matrix for the span 
			$$ D_i^{p} = \frac{p+1}{t_{i+p+1} - t_i} N_i^p $$
			
		where $N_i^p$ are the default splines of V space
	 	
	 	Require:
	 	--------
	 		@V is a FemTensor V= (V1, V2) with @V1, @V2 are FemSpace
	
		Returns:
		--------
			@M is a StencilMatrix - the mass matrix for D_i^p span	
	"""
	
	# ... sizes
	[s1, s2] = V.vector_space.starts
	[e1, e2] = V.vector_space.ends
	[p1, p2] = V.vector_space.pads
	# ...

	# ... seetings
	[k1, k2] = [W.quad_order for W in V.spaces]
	[spans_1, spans_2] = [W.spans for W in V.spaces]
	[basis_1, basis_2] = [W.quad_basis for W in V.spaces]
	[weights_1, weights_2] = [W.quad_weights for W in V.spaces]
	[points_1, points_2] = [W.quad_points for W in V.spaces]
	[V_1, V_2] = [W for W in V.spaces]

	# ... data structure
	M = StencilMatrix(V.vector_space, V.vector_space)
	# ...

	# ... build matrices
	for ie1 in range(s1, e1+1-p1):
		for ie2 in range(s2, e2+1-p2):
			i_span_1 = spans_1[ie1]
			i_span_2 = spans_2[ie2]
			for il_1 in range(0, p1+1):
				for jl_1 in range(0, p1+1):
					for il_2 in range(0, p2+1):
						for jl_2 in range(0, p2+1):
							i1 = i_span_1 - p1  - 1 + il_1
							j1 = i_span_1 - p1  - 1 + jl_1

							i2 = i_span_2 - p2  - 1 + il_2
							j2 = i_span_2 - p2  - 1 + jl_2

							v_s = 0.0
							for g1 in range(0, k1):
								for g2 in range(0, k2):
									di1_0 = coefficients_linear_combinaison(V_1, i1)*basis_1[ie1, il_1, 0, g1]
									di2_0 = coefficients_linear_combinaison(V_2, i2)*basis_2[ie1, il_2, 0, g1]
									bi_0 = di1_0 * di2_0
									
									d_j1 = basis_1[ie1, jl_1, 0, g1] * coefficients_linear_combinaison(V_1, j1)
									d_j2 = basis_2[ie2, jl_2, 0, g2] * coefficients_linear_combinaison (V_2, j2)
									bj_0 = d_j1 * d_j2
									
									wvol = weights_1[ie1, g1] * weights_2[ie2, g2]
									
									v_s = v_s + bi_0 * bj_0 * wvol
							M[i1, i2, j1 - i1, j2 - i2]  += v_s
	# ...
	return M
# ...
