__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-29"

# ...
from spl.core.bsplines import collocation_matrix

import numpy as np

from mat1d.discretdev import *
from space.spaces import *
from proj.pigrad import *
from proj.pil2 import *
from space.spaces import *

from spl.core.bsplines import collocation_matrix
# ..


if __name__ == "__main__":

	p = 2
	ne = 7
	grid = np.linspace(0., 1., ne+1)
	
	print('degree : ', p)
	print('number of points on the grid : ', ne)	
	
	V, W = creationSpaces(3, grid)
	V.init_collocation()
	
	
	print(V._bmat)
	print(V._bmat.shape)
	f = lambda x, t: np.cos(x)*np.cos(t)
	fprime = lambda x, t: - np.sin(x)*np.cos(t)
	
	
	col = collocation_matrix(V.knots, p, grid, False)
	print(col)
	print(col.shape)
	print('\n')
	
	
	pi0 = PiGrad (V, f, coef=False, dirichlet=True)
	pi1 = PiL2 (W, fprime, coef=True, dirichlet=True)
	
	G = DiscretDev(V, W)
	
	tmp = G.dot(pi0.proj())
	print(tmp)
	
	print(pi1.proj())
	
	#assert(np.allclose(tmp, rot.todense().transpose(), 1e-10))	
	#print(" -> SUCCESS")
# ...
	
	
	
	
# ...


