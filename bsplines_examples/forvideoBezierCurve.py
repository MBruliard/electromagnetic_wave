"""
	@author: Margaux BRULIARD
	@date: 18.06.2018
	@purpose: Example of a Bézier curve


	We know an analytical formula for a Bezier curve

	Order 1:  B(t) = (1-y)P0 + tP1
	
	Order 2:  B(t) =(1-t)²P0 + 2t(1-t)P1 + t²P2
"""


######################## MODULES ###########################

import numpy as np
from matplotlib import pyplot as plt



######################## USERS FUNCTION ##############################
def bezierOrder2 (t, x0, y0, x1, y1, x2, y2):
	
	xval = (1-t)**2*x0 + 2*t*(1-t)*x1 + t**2*x2
	yval = (1-t)**2*y0 + 2*t*(1-t)*y1 + t**2*y2


	return xval, yval

    
    
#==============================================================================
# MAIN FUNCTION
#==============================================================================
if __name__ == "__main__":

	x_min = -2.
	x_max = 2.
	p = 2;

	## We building control points
	x = np.zeros(p+1)
	y = np.zeros(p+1)

	x[0] = 0.
	x[1] = 0.3
	x[2] = -0.1

	y[0]  =0.
	y[1] = -0.3
	y[2] = 0.2

	for i in range (0, 100):
		y[1] = y[1] + 0.005
		x[1] = x[1] - 0.001
		## Bezier order 2 formula
		t = np.linspace(0,1,10)
		xval, yval = bezierOrder2 (t, x[0], y[0], x[1], y[1], x[2], y[2])


		plt.scatter(x, y, label='control points', color='r')
		plt.plot(xval, yval, label='Bezier curve', color='b')
		plt.title ('Bezier quadratic curve')
		plt.legend()
		plt.ylim(-0.35,  0.35)
		plt.xlim(-0.15, 0.35)
		plt.savefig('Images/video/beziercurve_it'+str(i)+'.png')
		print("The figure has been saved in 'Images/video/beziercurve_it"+str(i)+".png'")
		plt.clf()





