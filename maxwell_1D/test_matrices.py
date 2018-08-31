__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-21-08"

import numpy as np
from mat1d.mass import *
from mat1d.rotational import *
from mat1d.discretdev import *
from space.spaces import *
from userfunc.matrices_maxwell import assembly_rotational_matrix, assembly_derivate_discrete_matrix, assembly_mass_matrix_of_Dfamilly




if __name__ == "__main__":

	p = 3
	ne = 5
	grid = np.linspace(0., 1., ne)
	
	print('degree : ', p)
	print('number of points on the grid : ', ne)	
	
	V, W = creationSpaces(3, np.linspace(0., 1., ne))

	mass = MassMatrix (V)
	massD = MassMatrix(W, coef=True)
	rot = Rotational(V, W)
	dis = DiscretDev(V, W)

	tmp = np.dot(massD.todense(), dis.todense())
	assert(np.allclose(tmp, rot.todense().transpose(), 1e-10))	
	print(" -> SUCCESS")



