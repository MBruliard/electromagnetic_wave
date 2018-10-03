
from mat2d.rotational import *
from mat2d.mass import *
from mat2d.discretederivation import *
from space.spaces import *
from field import *


import numpy as np
from math import pi






if __name__ == "__main__":

	# Input data: degree, number of elements
	p1  = 2 ; p2  = 2
	ne1 = 3 ; ne2 = 3

	# Input analytical functions
	k1 = 1.
	k2 = 1.
	sigma1 = 0.
	sigma2 = 0.
	omega = pi
	
	
	analyticalh = lambda x, y, t: np.cos(k1 *x + sigma1)* np.sin(k2*y +sigma2)* np.cos(omega*t)
	analyicalex = lambda x,y, t: - k2/omega * np.cos(k1*x + sigma1)*np.sin(k2*y + sigma2) * np.sin(omega*t)
	analyticaley = lambda x, y, t: k2/omega * np.cos(k1*x + sigma1) * np.cos(k2*y + sigma2) * np.sin(omega*t)
	
	
	# Creating spaces
	V, Wdiv1, Wdiv2, W = creation_spaces_2d(p1, p2, ne1, ne2)

	# Building matrices
	mass  = Mass(V)
	mass = mass.todense()
	
	G = DiscreteDev (V)
	R = Rotational(V, Wdiv1, Wdiv2)

	valB = valuesfromgreville(V, analyticalh)
	valEx = valuesfromgreville(Wdiv1, analyicalex)
	valEy = valuesfromgreville(Wdiv2, analyticaley)

	startB = Field1d(valB)
	startE = Field2d(valEx, valEy)
	
	print("Size of vectors: ")
	print(startB.size())
	print(startE.dirx().shape)
	print(startE.diry().shape)
		
	print("\nTest of Rotational: ")
	print('shape R : ', R.Rs[0].shape, R.Rs[1].shape )
	test_calculR = R.dot(startE)
	print('shape of test_calculR: ', test_calculR.shape)
	
	print("\nTest of DiscretDev: ")
	print("shape G : ", G._dis[0].shape, G._dis[1].shape)
	test_calculG = G.dot(startB)
	print('shape of test_calculG: ', test_calculG.shape)
#	
#	tmp = massD.tocsc().dot()
#	assert(np.allclose(tmp, rot.todense().transpose(), 1e-10))	
#	
#	
#	
	
