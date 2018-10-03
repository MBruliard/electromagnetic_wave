"""
	@author: Margaux BRULIARD
	@date: 19.06.2018
	@purpose: A B-spline curve interpolation (basic example)
	
	
	Program:
		* knots, degree and control points are data
		* we compute the B-Spline family associated to the data
		* we compute the B-Spline curve
		
"""


from bsplines import Bspline
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import splev



#==============================================================================
# MAIN FUNCTION
#==============================================================================
print ("\nProgram Execution .... 'exampleBsplineCurve'")
choice_fig = input ("Do you want to save the graph (1) or do you prefer to plot it (2)  ?")
choice_fig = int(choice_fig)


# knots
knotlist = [0, 0, 1, 2, 3, 4, 5, 6, 6]
m = len(knotlist) #number knots
print('The list of knots: ', knotlist)
print("knots number: ", m)

 
x = [0, 0.5, 2.1, 3.4, 4  , 5.5, 6 ] #, 0.1, 0.2]
y = [0, 0.8, 0.5, 0.9, 2.3, 1.4, 0.6]
n = len(x)


# B-Spline curve degree
p = m - n 
print("B-spline curve degree: ", p)

if (not (n==len(x))):
	print("\nERROR: Controls points are false\n")
else:
	t = np.linspace (min(knotlist), max(knotlist), 20) 


	bsp = Bspline(knotlist, p); #Bspline family generating



	Sx = []
	Sy = []
	for tval in t:
		
		sum_x = 0
		sum_y = 0
		
		for i in range (0, n):
			bspi = bsp(tval, i=i)
			sum_x = sum_x + bspi*x[i]
			sum_y = sum_y + bspi*y[i]
		
		Sx.append(sum_x)
		Sy.append(sum_y)



	## Graphics 
	plt.scatter(x, y, label='control points', color='r')
	plt.plot(Sx, Sy, label='B-spline curve', color='b')
	plt.title ("B-splines curve interpolation based on control points")
	plt.legend()

	if choice_fig == 1 :
		plt.savefig('rapport/img/BsplineCurveGenerated.png');
		print("The figure has been save in 'rapport/img/BsplineCurveGenerated.png'")
	else:
		plt.show()

print("\n")






