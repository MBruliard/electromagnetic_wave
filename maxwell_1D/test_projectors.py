


from spl.core.bsplines import collocation_matrix

import numpy as np
from mat1d.mass import *
from mat1d.rotational import *
from mat1d.discretdev import *
from space.spaces import *

from proj.pigrad import InterpolationMatrix
from proj.pil2 import *

from space.spaces import *
from userfunc.matrices_maxwell import assembly_rotational_matrix, assembly_derivate_discrete_matrix, assembly_mass_matrix_of_Dfamilly




if __name__ == "__main__":

	p = 2
	ne = 7
	grid = np.linspace(0., 1., ne+1)
	
	print('degree : ', p)
	print('number of points on the grid : ', ne)	
	
	V, W = creationSpaces(3, grid)
	V.init_collocation()
	
#	f = lambda x, t: np.cos(x)*np.cos(t) 
#	
#	inter = InterpolationMatrix(V)
	knots = V.knots
	
	print('grid: ', grid)
	print('knots; ', knots)
	
	x =[0.] + list(np.linspace(0., 1., ne+1)) + [1.]
	m = collocation_matrix(knots, p, V.greville, False)
	print(m.shape)
	print(m)

#	histo = Histopolation(W, coeff=True)	
#	projL2 = PiL2(W, f, coeff=True )



