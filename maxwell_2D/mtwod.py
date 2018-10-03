__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-10"
__purpose__ = "Functions for maxwell 2D scheme "

from spl.fem.splines import SplineSpace
from spl.linalg.stencil import StencilVector, StencilMatrix
from spl.fem.tensor  import TensorFemSpace
import numpy as np


# --------------------------- ADDITIONNAL FUNCTIONS ------------------------ #
def coefficients_linear_combinaison(V, spline_index):
	"""
		Compute the alpha coefficient of the Schoenberg space as:
		
		$$ D_i^{p-1} = alpha_i * N_i^{p-1} $$
		
		with @V the schoenberg space of degree p-1 and @spline_index = i in this example 
	"""
	
	p = V.degree
	knots = V.knots

	if knots[spline_index + (p+1)] - knots[spline_index] == 0:
		return 0
	else:
		return (p+1)/(knots[spline_index + (p+1)] - knots[spline_index])


# ------------------- ASSEMBY MATRICES ------------------------------------ #
def assembly_mass_V(V):

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
    
#    print(V.spaces[0].quad_order)
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
                                    bi_0 = basis_1[ie1, il_1, 0, g1] * basis_2[ie1, il_2, 0, g2]
                                    bj_0 = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 0, g2]
#                                    bj_x = basis_1[ie1, jl_1, 1, g1] * basis_2[ie2, jl_2, 0, g2]
#                                    bj_y = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 1, g2]

                                    wvol = weights_1[ie1, g1] * weights_2[ie2, g2]

#                                    v_s += (bi_x * bj_x + bi_y * bj_y) * wvol
                                    v_s = v_s + bi_0 * bj_0 * wvol
                            M[i1, i2, j1 - i1, j2 - i2]  += v_s
    # ...
    return M
# ...


# ...
def assembly_mass_matrix_D(V):

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
    
#    print(V.spaces[0].quad_order)
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



def assembly_G_matrix (size_matrix):
	return 1










def assembly_R1_matrix (V, Wdiv, sizeInTuple):
	"""
		Require:
		--------
			@V is a FemSpace represents H1 -> V = (V1, V2)
			@Wdiv is a FemSpace represents H1 curl space -> Wdiv = (V1, W2)
			@sizeInTuple is a tuple (i1, i2, j1, j2) represents the size of the matrix R
		
		Returns:
		--------
			@R1 a matrix of size (i1, i2, j1, j2) 
	"""

	R1 = np.ndarray(sizeInTuple)

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
	[kw1, kw2] = [W.quad_order for W in Wdiv.spaces]
	[spans_w1, spans_w2] = [W.spans for W in Wdiv.spaces]
	[basis_w1, basis_w2] = [W.quad_basis for W in Wdiv.spaces]
	#    [weights_w1, weights_w2] = [W.quad_weights for W in Wdiv.spaces]
	#    [points_1, points_2] = [W.quad_points for W in V.spaces]
	[Wdiv_1, Wdiv_2] = [W for W in Wdiv.spaces]

	print("k1 et k2: ", k1, k2)
	print("kw1 et kw2: ", kw1, kw2)
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
									
									vm = vm + (bj1_0 * dj2_0 * bi1_0 * bi2_x)
								# ...
							# ...
							R1[i1, i2, j1 - i1, j2 - i2] = R1[i1, i2, j1 - i1, j2 - i2] + vm					
						# ...
					# ...
				#...
			# ...
		# ...
	# ...
	return R1




# .....
def assembly_R2_matrix (Vspace, Wspace, sizeInTuple ):
	"""
		Require:
		--------
			@V is a FemSpace represents H1 -> V = (V1, V2)
			@Wdiv is a FemSpace represents H1 curl space -> Wdiv = (W1, V2)
			@sizeInTuple is a tuple (i1, i2, j1, j2) represents the size of the matrix R2
		
		Returns:
		--------
			@R2 a matrix of size (i1, i2, j1, j2) 
	"""
	
	R2 = np.ndarray(sizeInTuple)
		
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
			
			for il_1 in range(0, pv1+1):
				i1 = i_spans_1 - pv1 + il_1
			
				for jl_1 in range(0, pw1 + 1):
					j1 = j_spans_1 - pw1 + jl_1
					
					for il_2 in range(0, pv2+1):
						i2 = i_spans_2 - pv2 + il_2
						
						for jl_2 in range (0, pw2):
							j2 = i_spans_2 - pw2 + jl_2
							
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
							R2[i1, i2, j1 - i1, j2 - i2] = R2[i1, i2, j1 - i1, j2 - i2] + vm					
						# ...
					# ...
				#...
			# ...
		# ...
	# ...
	return R2


def assembly_R_matrix (V, Wdiv1, Wdiv2):


	R1 = assembly_R1_matrix(V, Wdiv, size1)
	R2 = assembly_R2_matrix(V, Wdiv, size2)
	
	R = [R1, R2]

	 




