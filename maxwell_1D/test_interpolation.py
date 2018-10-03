__author__ = "Margaux BRULIARD"
__date__ = "2018-08-02"
__partners__ = "Ahmed RATNANI", "NMPP"
__purpose__ = "test of the assembly of matrices_maxwell.py"

from math import pi

from userfunc.matrices_maxwell import *
from userfunc.projector import *
from userfunc.norm import *
from userfunc.moned import *
from userfunc.energy import *

if __name__ == "__main__":

	# -- Data	
	number_of_elements = 5
	grid =  np.linspace(0., 1., number_of_elements)
	degree = 2

	
	omega = 1.
	k = - omega
	
	mu0 = 1.
	epsilon0 = 1.
	
	analytical_e = lambda x, t: np.sin(2*pi*k*x)*np.cos(omega*t)
	analytical_b = lambda x, t: np.cos(2*pi*k*x)*np.sin(omega*t)
	
	V, W = creation_spaces(degree, grid)
	
	[s] = V.vector_space.starts
	[e] = V.vector_space.ends
	[p] = V.vector_space.pads
	
	spans = V.spans
	basis = V.quad_basis
	
	print(spans)
	
	for ie1 in range(s, e-p+1):
		i_span = spans[ie1]
		
		
		for il in range (0, p):
			i = i_span - p + il
			
			print('ie1='+str(ie1)+' spline '+str(i) + ' x0 = '+str(basis[ie1, il, 0, 0])+ ' xn='+str(basis[ie1, il, 0, -1]))	
	
	
	
	
	
	
	
	
	
	 
