
__author__ = "Margaux BRULIARD"
__date__ = "2018-07-12"
__partners__ = "Ahmed RATNANI", "NMPP"
__purpose__ = "mass matrix maxwell 1D"


# -- spl packages importing
from spl.fem.splines import SplineSpace 


def creationSpaces (degree, grid):
	"""
        Create the schoenberg spaces of degree @degree (and @degree-1)
        
        Parameters
        ----------
        @degree: int
			degree of the first schoenberg space         	

        @grid : list of float or numpy.ndarray
           	breakpoints of the domain

        Returns
        -------
        V : Class SplineSpace
            1D finite element space of schoenberg space of degree @degree
        
        W: Class SplineSpace
        	1D finite element space of schoenberg space of degree (@degree - 1)

        """


	V = SplineSpace(degree, grid=grid)
	W = SplineSpace(degree-1, grid=grid)

	V.init_fem(quad_order=degree+1)
	W.init_fem(quad_order=degree+1)

	return V, W

# ...


def coeffDspan(V, splineIndex):
	"""
		Compute the alpha coefficient of the Schoenberg space as:
		
		$$ D_i^{p-1} = alpha_i * N_i^{p-1} $$
		
		with @V the schoenberg space of degree p-1 and @spline_index = i in this example 
	"""
	
	p = V.degree
	knots = V.knots

	if knots[splineIndex + (p+1)] - knots[splineIndex] == 0:
		return 0
	else:
		return (p+1)/(knots[splineIndex + (p+1)] - knots[splineIndex])

# ...


def getNumberOfSplines (V):
	"""
		Require:
		--------
			V:  SplineSpace
			
		Returns:
		--------
			nbS: Integer number of Splines in V
	"""
	
	knots = V.knots
	
	nbS = len(knots) - V.degree -1
	
	return nbS







