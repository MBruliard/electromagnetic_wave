__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-22"
__purpose__ = "assembly rotational matrix for Maxwell 2D equation"

from spl.fem.splines import SplineSpace
from spl.linalg.stencil import StencilVector, StencilMatrix
from spl.fem.tensor  import TensorFemSpace

import numpy as np
from scipy.sparse import csc_matrix, csr_matrix

from space.spaces import coefficients_linear_combinaison, getNumberSplineInSpace
from util.mattransform import *
from field import *

class Rotational:
	"""
		Compute the rotational matrix for the solver of Maxwell 2D equation
		
		Entries:
		--------
			V:  space which represents H^1 space
			
			Wdiv1 and Wdiv 2: reprensent (H^1, curl) space
			
		Attributs:
		----------
			Rs: list of 2 sparse matrix 2D -> Rs[0] = R1 and Rs[1] = R2
	"""
	
	def __init__ (self, V, Wdiv1, Wdiv2):
		
		# --  compute size of R1 and R2
		[V1, V2] = [W for W in V.spaces]
		[Wdiv1_1, Wdiv1_2] = [W for W in Wdiv1.spaces]
		[Wdiv2_1, Wdiv2_2] = [W for W in Wdiv2.spaces]
		
		nbs_v1 = getNumberSplineInSpace(V1)
		nbs_v2 = getNumberSplineInSpace(V2)	
		nbs_w11 = getNumberSplineInSpace(Wdiv1_1)
		nbs_w12 = getNumberSplineInSpace(Wdiv1_2)
		nbs_w21 = getNumberSplineInSpace(Wdiv2_1)
		nbs_w22 = getNumberSplineInSpace(Wdiv2_2)
		
		print('shape V: ', nbs_v1, nbs_v2)
		print('shape Wdiv1: ', nbs_w11, nbs_w12)
		print('shape Wdiv2: ', nbs_w21, nbs_w22)
		# -- compute R1 & R2
		R1 = kernelR1(V, Wdiv1, (nbs_v1, nbs_v2, nbs_w11, nbs_w12))
		R2 = kernelR2(V, Wdiv2, (nbs_v1, nbs_v2, nbs_w21, nbs_w22))
		
		
		self.Rs = [R1, R2]
		self.size = (R1.shape, R2.shape)
	# ...
	
	
	def dot (self, xs):
		"""
			Compute R*xs
			
			Entries:
			--------
				xs: Field2d
				
			Output:
			-------
				y: Field1d
		"""
		
		x1 = xs.dirx()
		x2 = xs.diry()
		

##		y1 = nonameyet(self.Rs[0], x1)
##		y2 = nonameyet(self.Rs[1], x2)
#		
#		return Field1d(y1+y2)

		print('test shpae of R1, x1 et x2: ', self.Rs[0].shape, x1.shape, x2.shape)
		print('shape of R2: ', self.Rs[1].shape)
		y1 = self.Rs[0].dot(x1)
		y2 = self.Rs[1].dot(x2)
		print("shape y1 in field: ", y1.shape)
		print("shape y2 in field: ", y2.shape)
		return Field1d(y1 + y2)

			
	# ... 
	
	def printf(self):
		
		print([self.Rs[0], self.Rs[1]])
	# ...
	
	
#	def todense(self):
#		
#		self.Rs[0] = self.Rs[0].todense()
#		self.Rs[1] = self.Rs[1].todense()	
#	# ..
	



def kernelR1 (V, Wdiv, sizeInTuple):
	"""
		Require:
		--------
			@V is a FemSpace represents H1 -> V = (V1, V2)
			@Wdiv is a FemSpace represents H1 curl space -> Wdiv = (V1, W2)
			@sizeInTuple is a tuple (i1, i2, j1, j2) represents the size of the matrix R
		
		Returns:
		--------
			@R1 a sparse matrix of size (i1*j1, i2*j2) 
	"""

	R1 = np.zeros(sizeInTuple)

	# ... sizes
	[sv1, sv2] = V.vector_space.starts
	[ev1, ev2] = V.vector_space.ends
	[pv1, pv2] = V.vector_space.pads
	# ...

	# ... seetings
	[k1, k2] = [W.quad_order for W in V.spaces]
	[spans_v1, spans_v2] = [W.spans for W in V.spaces]
	[basis_v1, basis_v2] = [W.quad_basis for W in V.spaces]
	[weights_v1, weights_v2] = [W.quad_weights for W in V.spaces]
	[V_1, V_2] = [W for W in V.spaces]

	# ... sizes
	[sw1, sw2] = Wdiv.vector_space.starts
	[ew1, ew2] = Wdiv.vector_space.ends
	[pw1, pw2] = Wdiv.vector_space.pads
	# ...

	# ... seetings
	[kw1, kw2] = [W.quad_order for W in Wdiv.spaces]
	[spans_w1, spans_w2] = [W.spans for W in Wdiv.spaces]
	[basis_w1, basis_w2] = [W.quad_basis for W in Wdiv.spaces]
	[Wdiv_1, Wdiv_2] = [W for W in Wdiv.spaces]


	# .... start
	for ie1 in range(sv1, ev1 - pv1 + 1 ):
		spans_1 = spans_v1[ie1]
		
		for ie2 in range(sv2, ev2 - pv2 + 1):
			i_spans_2 = spans_v2[ie2]
			j_spans_2 = spans_w2[ie2]
			
			for il_1 in range(0, pv1+1):
				i1 = spans_1 - pv1 + il_1
			
				for jl_1 in range(0, pw1 + 1):
					j1 = spans_1 - pw1 + jl_1
					
					for il_2 in range(0, pv2+1):
						i2 = i_spans_2 - pv2 + il_2
						
						for jl_2 in range (0, pw2):
							j2 = j_spans_2 - pw2 + jl_2
							
							vm = 0.0
							for g1 in range (0, k1):
								for g2 in range (0, k2):
									bj1_0 = basis_w1[ie1, jl_1, 0, g1]
									dj2_0 = basis_w2[ie2, jl_2, 0, g2]*coefficients_linear_combinaison(Wdiv_2, j2)
									
									bi1_0 = basis_v1[ie1, il_1, 0, g1]
									bi2_x = basis_v2[ie2, il_2, 1, g2]
									
									wvol = weights_v1[ie1, g1] * weights_v2[ie2, g2]
									
									vm = vm + (bj1_0 * dj2_0 * bi1_0 * bi2_x)* wvol
								# ...
							# ...
							R1[i1, i2, j1, j2 ] = R1[i1, i2, j1, j2 ] + vm					
						# ...
					# ...
				#...
			# ...
		# ...
	# ...
	return R1
#	return csr_matrix(fourToTwoD(R1))
# ...

def kernelR2 (V, Wdiv, sizeInTuple ):
	"""
		Require:
		--------
			@V is a FemSpace represents H1 -> V = (V1, V2)
			@Wdiv is a FemSpace represents H1 curl space -> Wdiv = (W1, V2)
			@sizeInTuple is a tuple (i1, i2, j1, j2) represents the size of the matrix R2
		
		Returns:
		--------
			@R2 a sparse matrix of size (i1*j1, i2*j2) 
	"""
	# TODO revoir la construction de cette matrice !!!!
	R2 = np.zeros(sizeInTuple)
		
	# ... sizes
	[sv1, sv2] = V.vector_space.starts
	[ev1, ev2] = V.vector_space.ends
	[pv1, pv2] = V.vector_space.pads
	# ...

	# ... seetings
	[k1, k2] = [W.quad_order for W in V.spaces]
	[spans_v1, spans_v2] = [W.spans for W in V.spaces]
	[basis_v1, basis_v2] = [W.quad_basis for W in V.spaces]
	[weights_v1, weights_v2] = [W.quad_weights for W in V.spaces]
	#    [points_1, points_2] = [W.quad_points for W in V.spaces]
	[V_1, V_2] = [W for W in V.spaces]

	# ... sizes
	[sw1, sw2] = Wdiv.vector_space.starts
	[ew1, ew2] = Wdiv.vector_space.ends
	[pw1, pw2] = Wdiv.vector_space.pads
	# ...

	# ... seetings
	[spans_w1, spans_w2] = [W.spans for W in Wdiv.spaces]
	[basis_w1, basis_w2] = [W.quad_basis for W in Wdiv.spaces]
	#    [weights_w1, weights_w2] = [W.quad_weights for W in Wdiv.spaces]
	#    [points_1, points_2] = [W.quad_points for W in V.spaces]
	[Wdiv_1, Wdiv_2] = [W for W in Wdiv.spaces]

	
	# .... start
	for ie1 in range(sv1, ev1 - pv1 + 1 ):
		i_spans_1 = spans_v1[ie1]
		j_spans_1 = spans_w1[ie1]
		for ie2 in range(sv2, ev2 - pv2 + 1):
			i_spans_2 = spans_v2[ie2]
			j_spans_2 = spans_w2[ie2]
			
			for il_1 in range(0, pv1+1):
				i1 = i_spans_1 - pv1 + il_1
			
				for jl_1 in range(0, pw1 + 1):
					j1 = j_spans_1 - pw1 + jl_1
					
					for il_2 in range(0, pv2+1):
						i2 = i_spans_2 - pv2 + il_2
						
						for jl_2 in range (0, pw2+1):
							j2 = j_spans_2 - pw2 + jl_2
							
							vm = 0.0
							for g1 in range (0, k1):
								for g2 in range (0, k2):
									dj1_0 = basis_w1[ie1, jl_1, 0, g1]* coefficients_linear_combinaison(Wdiv_1, j1)
									bj2_0 = basis_w2[ie2, jl_2, 0, g2]
									
									bi1_x = basis_v1[ie1, il_1, 1, g1]
									bi2_0 = basis_v2[ie2, il_2, 0, g2]
									
									wvol = weights_v1[ie1, g1] * weights_v2[ie2, g2]
									
									vm = vm + (- dj1_0 * bj2_0 * bi1_x * bi2_0)
								# ...
							# ...
							R2[i1, i2, j1, j2] = R2[i1, i2, j1, j2] + vm					
						# ...
					# ...
				#...
			# ...
		# ...
	# ...
	return R2
#	return csr_matrix(fourToTwoD(R2))


