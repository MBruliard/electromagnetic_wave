__author__ = "Margaux BRULIARD"
__date__ = "20.06.2018"
__partners__ = "Ahmed RATNANI", "NMPP from IPP Munchen"
__purpose__ = "Solver Poisson Equation 1D with homogenous dirichlet boundary conditions using SPL Python package"



# --- Importing modules --- #
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve


from spl.fem.splines import SplineSpace 
from spl.utilities.quadratures import gauss_legendre
from spl.linalg.stencil        import StencilVector, StencilMatrix

from matplotlib import pyplot as plt

class StiffnessMatrix:
	"""
		define the stiffness Matrix of a SplineSpace 1D
		
		Attribut:
		---------
			_stiff : StencilMatrix
	"""

	def __init__ (self, V):
		"""
			Entry:
			------
				V: SplineSpace used
		"""
	# ... sizes
		[s1] = V.vector_space.starts
		[e1] = V.vector_space.ends
		[p1] = V.vector_space.pads

		k1        = V.quad_order
		spans_1   = V.spans
		basis_1   = V.quad_basis
		weights_1 = V.quad_weights
		# ...

		# ... data structure
		stiffness = StencilMatrix( V.vector_space, V.vector_space )
		# ...

		# ... build matrices
		for ie1 in range(s1, e1+1-p1):
		    i_span_1 = spans_1[ie1]
		    for il_1 in range(0, p1+1):
		        for jl_1 in range(0, p1+1):
		            i1 = i_span_1 - p1 + il_1
		            j1 = i_span_1 - p1 + jl_1

		            v_s = 0.0
		            for g1 in range(0, k1):
		                bi_x = basis_1[ie1, il_1, 1, g1]

		                bj_x = basis_1[ie1, jl_1, 1, g1]

		                wvol = weights_1[ie1, g1]

		                v_s += (bi_x * bj_x) * wvol

		            stiffness[i1, j1 - i1]  += v_s
		# ...

		# ...
		
		self._stiff = stiffness	
	# ...
	
	def todense(self):
		return self._stiff.tocoo().todense()
#		return self._stiff
	# ...
	
	def tocsc(self):
		return csc_matrix( self._stiff.tocoo() )
	# ...
	
	def shape(self):
		return self.todense().shape
	# ...
# ...


def assemblyRhs (V, func):
	"""
		Build the second member of our system
		
		Entries:
		--------
			V: SplineSpace
			
			func: second member function func(x)
			
		Output:
		-------
			rhs: np.ndarray 2D (n, 1) 
	"""
	# ... sizes
	[s1] = V.vector_space.starts
	[e1] = V.vector_space.ends
	[p1] = V.vector_space.pads

	k1        = V.quad_order
	spans_1   = V.spans
	basis_1   = V.quad_basis
	points_1  = V.quad_points
	weights_1 = V.quad_weights
	# ...

	# ... data structure
#	rhs = StencilVector( V.vector_space )
	rhs = np.zeros((V.nbasis, 1))
	# ...

	# ... build rhs
	for ie1 in range(s1, e1+1-p1):
		i_span_1 = spans_1[ie1]
		for il_1 in range(0, p1+1):
			i1 = i_span_1 - p1 + il_1

			v = 0.0
			for g1 in range(0, k1):
				bi_0 = basis_1[ie1, il_1, 0, g1]
				bi_x = basis_1[ie1, il_1, 1, g1]

				x1    =  points_1[ie1, g1]
				wvol  = weights_1[ie1, g1]

				#                v += bi_0 * np.sin( 2*np.pi*x1 )*(2*np.pi)**2 * wvol
				v += bi_0 * func(x1) * wvol

			rhs[i1] += v
	# ...
	return rhs
# ...

def computeAnalayticalSolution(abslist, func):

	y = np.ndarray(abslist.shape)

	for i in range(len(abslist)):
		y[i] = func (abslist[i])
		
	return y
# ...


if __name__ == "__main__":
	
	# defintiion of functions
	sec = lambda x: np.sin( 2*np.pi*x )*(2*np.pi)**2
	analytic = lambda x: np.sin(2*np.pi*x)
	
	# define dirichlet boundary conditions
	dirichletleft = 0.
	dirichletright = 0.
	
	# SplineSpace
	p = 3
	ne = 60
	grid = np.linspace(0., 1., ne+1)

	
	# Main part
	print("Programm running ... ")
	
	
	V = SplineSpace(p, grid=grid) ; V.init_fem()
	
	# building stiffness matrix in a sparse matrix
	stiff = StiffnessMatrix(V)
	stiff = stiff.todense()


	rhs = assemblyRhs(V, sec)


	# homogenous dirichlet boundary conditions
	stiff = stiff[1:-1, 1:-1]
	rhs = rhs[1:-1]

	u = np.linalg.solve(stiff, rhs)
	u = list(u)
	u.insert(0, dirichletleft)
	u.append(dirichletright)


	greville = V.greville
	uex = computeAnalayticalSolution(greville, analytic)
	
	# ploting results
	plt.plot(greville, u, label='Approximate solution')
	plt.plot (greville, uex, label='analytic solution')
	plt.legend()
	plt.title('Example of Poisson 1D solver\np='+str(p)+' ne='+str(ne))
	plt.show()
	
	
	# saving graph
	savefig = input ("Do you want to save this graphic ? [O/n] : ")
	if (savefig == 'O'):
		plt.savefig('Images/poisson1dspl_ne'+str(ne)+'_p'+str(p)+'.png')
		print("Graphic saved in 'Images/poisson1dspl_ne'+str(ne)+'_p'+str(p)+'.png'")
	
# ...
