#from mtwod import *
from mat2d.mass import *
from space.spaces import *
from mat2d.rotational import *
from kronecker import *
from mattransform import *
import numpy as np


# ....
if __name__ == "__main__":

		
	# Input data: degree, number of elements
	p1  = 3  ; p2  = 3
	ne1 = 5 ; ne2 = 5

	V, Wdiv1, Wdiv2, W = creation_spaces_2d(p1, p2, ne1, ne2)


	# -- test de la fonction kronecker
	A= np.random.rand(2, 2)
	B = np.random.rand(2, 2)
	
	K = kronecker(A, B, twodim='false')
	
	Knew = mat_4d_to_2d(K)
	print("\nTest mattransform: ")
	print(Knew)


	
# ....



## ... numbers of elements and degres
#p1  = 2 ; p2  = 2
#ne1 = 3 ; ne2 = 3
## ...

#print('> Grid   :: [{ne1},{ne2}]'.format(ne1=ne1, ne2=ne2))
#print('> Degree :: [{p1},{p2}]'.format(p1=p1, p2=p2))

#grid_1 = linspace(0., 1., ne1+1)
#grid_2 = linspace(0., 1., ne2+1)

#V1 = SplineSpace(p1, grid=grid_1) ; V1.init_fem()
#V2 = SplineSpace(p2, grid=grid_2) ; V2.init_fem()

#V = TensorFemSpace (V1, V2)
#M = assembly_v0(V)
#M = M.tocoo().todense()

#print(M)

