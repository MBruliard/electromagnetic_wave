__author__ = "Margaux BRULIARD"
__date__ = "20.06.2018"
__partners__ = "Ahmed RATNANI", "NMPP", "IPP"
__purpose__ = "Solver Poisson Equation 1D with homogenous dirichlet boundary conditions"
	

# importing modules
import numpy as np
from scipy.integrate import quadrature, quad
from matplotlib import pyplot as plt

from bsplines import Bspline, make_knots



def analyticSol (x):
	"""
		Analtic solution of our problem
	"""
	return x*(1-x);
# ...




def source (x):
	"""
		source term function
	"""
	return 2;
# ...


def legendreGauss (func, deg, a, b, ind, bsp, ind2=0):

	"""
		Compute the quadrature Gauss-Legendre method.
		
		Entries:
		--------
			func: Function we are studying
			
			deg: degree for Gauss Legendre quadrature method
			
			a, b : interval of integration
			
			ind: index of bspline
			
			bsp: list of bsplines
		
		Output:
		-------
			gauss: approximation of the integrale of func on [a, b] interval
	"""

	x, w = np.polynomial.legendre.leggauss(deg)
	t = 0.5*(x+1)*(b-a)+ a
	
	gauss = sum(w + func(t, bsp, ind, ind2))*( 0.5*(b-a))

	return gauss
# ...

	
def secondMember(x, bsp, i, j=0):
	"""
		Compute the second member of our system.
		
		Entries:
		---------
			x: variable to integrate
			
			bsp: B-splines family
			
			i: index of B-spline
		
		Output:
		-------
			y : values of the second member
	"""
	y = source(x)*bsp(x, i=i)
	return y
# ...

	
def bilinearForm (x, bsp, i, j):
	"""
		Function which represents the bilinar form of the discretize equation : \nabla N_j \nabla N_i $$
		
		Entries:
		--------
			x: float - variable to integrate
			
			bsp: Bspline family
			
			i, j: int - index of Bspline basis
			
		Output:
		------- 
			y: float - value of \nabla N_j(x) \nabla N_i(x)
	"""

	
	return bsp(x, i=i, n_deriv=1)*bsp(x, i=j, n_deriv=1);


def discretizeSecondMember (bsp, knotlist,p, nbquadrature):
	"""
		Discretize the second member of our equation to build the system
		
		Entries:
		--------
			bsp: Bspline Family
			
			knotlist: list of knots
			
			nbquad: quadrature degree for the gauss legendre method
		
		Output:
		-------
			F: dense matrix  
	"""
	F = np.zeros(bsp.N-2)
	for i in range(1, bsp.N-1): #i est l'indice des Ni
		
		
		# on regarde sur tous les knots pour chaque Ni
		for iknot in range(len(knotlist)-1):
			
			F[i-1] = F[i-1] + legendreGauss(secondMember, nbquadrature, knotlist[iknot], knotlist[iknot+1], i, bsp)
					
	return F;
# ...


def stiffnessMatrix (bsp, knotlist,p, nbquad):
	"""
		Compute the stiffness matrix
		
		Entries:
		--------
			bsp: Bsplines family
			
			knotlist: list of knots
			
			p: degree of B splines family
			
			nbquad: quadrature degree of the Gauss Legendre approximation
			
		Ouput:
		------
			S: dense matrix - stiffness matrix
	"""

	S = np.zeros((bsp.N-2, bsp.N-2))
	
	for line in range(1, bsp.N-1):
		
		for column in range(1, bsp.N-1):
		
			for iknot in range(len(knotlist)-1):
				S[line-1, column-1] = S[line-1, column-1] + legendreGauss(bilinearForm, nbquad, knotlist[iknot], knotlist[iknot+1], line, bsp, ind2=column)
	
	return  S;
	
# ---------------------------------- MAIN FUNCTION --------------------------------------- #
if __name__ == "__main__":

	# degree
	p = 2
	ne = 9
	knotlist = make_knots(ne, p)

	bsp = Bspline(knotlist, p);


	n = len(knotlist) - p # number of control points
	m = len(knotlist) - p - 1 # number of B-Splines

	# building matrices
	F = discretizeSecondMember (bsp, knotlist, p, 300)
	S = stiffnessMatrix (bsp, knotlist, p, 300)
	 

	a = knotlist[0]
	b = knotlist[-1]
	
	# analytic solution computing
	t = np.linspace(a, b, 20)
	uex = analyticSol (t)

	# solving 
	u = np.linalg.solve(S, F)

	# homogenous Dirichlet boundary conditions
	ufinal = np.insert(u, 0, 0)
	ufinal = np.append(ufinal, 0)


	# ploting 
	plt.plot(np.linspace(a, b, m), ufinal, label='solution u', color='b')
	plt.plot(t, uex, label='analytical solution', color='r')
	plt.legend()
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Approximate solution u based on B-splines\np='+str(p)+' ne='+str(ne))
	plt.show()

	savefig = input("Do you want to saved this graphic ? [O/n]: ")
	if (savefig == 'O'):
		plt.savefig('Images/poisson1d_ne'+str(ne)+'_p'+str(p)+'.png')
		print("The figure has been saved in 'Images/poisson1d_ne'+str(ne)+'_p'+str(p)+'.png'")



