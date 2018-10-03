__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-19"
__purpose__ = "user functions about spaces in Maxwell 2D"


from spl.fem.splines import SplineSpace
from spl.linalg.stencil import StencilVector, StencilMatrix
from spl.fem.tensor  import TensorFemSpace
import numpy as np

def coefficients_linear_combinaison(V, spline_index):
	"""
		Compute the alpha coefficient of the Schoenberg space as:
		
		$$ D_i^{p-1} = alpha_i * N_i^{p-1} $$
		
		with @V the schoenberg space of degree p-1 and @spline_index = i in this example 
		
		Require:
		--------
			@V is a FemSpace
			@spline_index is the index of the spline we are studying
		
		Returns:
		--------
			alpha_i coefficient presented before - float
			
	"""
	
	p = V.degree
	knots = V.knots

	if knots[spline_index + (p+1)] - knots[spline_index] == 0:
		return 0
	else:
		return (p+1)/(knots[spline_index + (p+1)] - knots[spline_index])
		


# -----------		
def getNumberSplineInSpace(V):
	"""
		Require:
		--------
			@V a FemSpace
			
		Returns:
		--------
			@nbS number of Splines in V space
	"""
	
	knots = V.knots
	
	nbS = len(knots) - V.degree -1
	
	return nbS
	
	
	
# -----------------
def creation_spaces_2d(p1, p2, ne1, ne2):
	"""
		Creation of the different spaces used in our solver:
		
		Require:
		--------
			@p1 and @p2 are degrees of direction 1 and 2
			@ne1 and @ne2 are number of elements in direction 1 and 2
		
		Returns:
		--------
			@V TensorFem which represents H^1 space
			@W TensorFem which represents L^2 space
			@Wdiv1 and @Wdiv2 represent (H¹, div) space
			
	"""
	grid_1 = np.linspace( 0., 1., num=ne1+1 )
	grid_2 = np.linspace( 0., 1., num=ne2+1 )

	# Create 1D finite element spaces and precompute quadrature data
	V1 = SplineSpace( p1, grid=grid_1); V1.init_fem(quad_order=p1+1)
	W1 = SplineSpace(p1-1, grid=grid_1); W1.init_fem(quad_order=p1+1)
	V2 = SplineSpace( p2, grid=grid_2 ); V2.init_fem(quad_order=p2+1)
	W2 = SplineSpace(p2-1, grid=grid_2); W2.init_fem(quad_order=p2+1)

	# Create 2D tensor product finite element space
	V = TensorFemSpace(V1, V2 )
	W = TensorFemSpace(W1, W2)

	Wdiv1 = TensorFemSpace(V1, W2)
	Wdiv2 = TensorFemSpace(W1, V2)
	
	return V, Wdiv1, Wdiv2, W







