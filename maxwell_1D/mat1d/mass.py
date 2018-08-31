__author__ = "Margaux BRULIARD"
__date__ = "2018-07-12"
__partners__ = "Ahmed RATNANI", "NMPP"
__purpose__ = "mass matrix maxwell 1D"


# -- spl packages importing
from spl.fem.splines import SplineSpace 
from spl.utilities.quadratures import gauss_legendre
from spl.linalg.stencil        import StencilVector, StencilMatrix

# ---
from space.spaces import coeffDspan

# -- 
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import inv


class MassMatrix:
	
	
	def __init__(self, V, coef=False):
		"""
			Assembly the mass matrix of the space V
			
			Entries:
			--------
				
				V: SplineSpace
				
				coef: boolean which choose (N_i^p) spans or (D_i^p) = (\frac{p+1}{t_{i+p+1} - t_i}N_i^p) spans
			
				dirichlet: boolean 
			Output:
			-------
				
				_mass : StencilMatrix - the mass matrix		
		"""	
		
		if (coef):
			self._mass = assemblyDSpanMassMatrix(V)
			
		else:
			self._mass = assemblyMassMatrix(V)
		
			
		self._size = self._mass.shape
	
	# ...
	
	def tocsc(self):
		return csc_matrix(self._mass.tocoo())
	# ...
	
	def dirichlet(self):
		return self.tocsc()[1:-1, 1:-1]
	
	def inv (self):
		"""
			Returns the inverse matrix in a sparse matrix
		"""
		masscsc = csc_matrix(self._mass.tocoo())
		return inv(masscsc)
	
	# ...
	
	def dot (self, mat):
		"""
			Compute the dot between the mass matrix and mat (vector or matrix)
			
			Return a sparse matrix (or vector)
		"""
		mass = self._mass.tocsc()
		
		return mass.dot(mat)
		
		
		
	# ...	
	
	def todense(self):
		"""
			Returns the dense matrix equivalent to _mass attribut
		"""
		return self._mass.tocoo().todense()
		
		
	def size(self):
		return self._size
# ...


def assemblyDSpanMassMatrix(V):
	"""
		Assembly the mass matrix of the space V
		
		Entries:
		--------
			
			V: SplineSpace
		
		Output:
		-------
			
			mass : StencilMatrix - the mass matrix
	"""
	# ... sizes
	[s1] = V.vector_space.starts
	[e1] = V.vector_space.ends
	[p1] = V.vector_space.pads

	k1        = V.quad_order
	spans_1   = V.spans
	basis_1   = V.quad_basis
	weights_1 = V.quad_weights
	knots_1 = V.knots
	# ...

	# ... data structure
	mass      = StencilMatrix( V.vector_space, V.vector_space )
	# ...

	# ... build matrices
	for ie1 in range(s1, e1+1-p1):
		i_span_1 = spans_1[ie1]
		for il_1 in range(0, p1+1):
		    for jl_1 in range(0, p1+1):
		        i1 = i_span_1 - p1 + il_1
		        j1 = i_span_1 - p1 + jl_1

		        # -- assembly D_j^p B-spline family: D_j^p = alpha * basis_w_0
		        alpha_j = coeffDspan(V, j1)               

		        # -- assembly D_i^p B-spline family: D_i^p = alpha * basis_w_0
		        alpha_i = coeffDspan(V, i1)

		        v_m = 0.0
		        for g1 in range(0, k1):
		        	bi_0 = basis_1[ie1, il_1, 0, g1]
		        	bj_0 = basis_1[ie1, jl_1, 0, g1]
		        	di_0 = bi_0 * alpha_i
		        	dj_0 = bj_0 * alpha_j
		        	wvol = weights_1[ie1, g1]
		        	v_m += di_0 * dj_0 * wvol
		            
		        mass[i1, j1 - i1] += v_m
	return mass
# ...


def assemblyMassMatrix(V):
	"""
		Assembly the mass matrix of the space V
			
			Entries:
			--------
				
				V: SplineSpace
			
			Output:
			-------
				
				mass : StencilMatrix - the mass matrix
	"""
	
	# ... sizes
	[s1] = V.vector_space.starts
	[e1] = V.vector_space.ends
	[p1] = V.vector_space.pads

	k1        = V.quad_order
	spans_1   = V.spans
	basis_1   = V.quad_basis
	weights_1 = V.quad_weights
	# ...

	# ... data structure
	mass      = StencilMatrix( V.vector_space, V.vector_space )
	# ...

	# ... build matrices
	for ie1 in range(s1, e1+1-p1):
		i_span_1 = spans_1[ie1]
		for il_1 in range(0, p1+1):
		    for jl_1 in range(0, p1+1):
		        i1 = i_span_1 - p1 + il_1
		        j1 = i_span_1 - p1 + jl_1

		        v_m = 0.0
		        for g1 in range(0, k1):
		            bi_0 = basis_1[ie1, il_1, 0, g1]

		            bj_0 = basis_1[ie1, jl_1, 0, g1]
		            
		            wvol = weights_1[ie1, g1]

		            v_m += bi_0 * bj_0 * wvol
		            
		        mass[i1, j1 - i1] += v_m
	# ...

	# ...
	return mass
# ...








