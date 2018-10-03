__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-12"
__aim__ = "Projectors used into DeRham sequence for maxwell 1D"


# --- Importing modules --- #
import numpy as np
from sys import exit
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

from mat1d.mass import *
from space.spaces import getNumberOfSplines
from field import Field

from spl.fem.splines import SplineSpace 
from spl.utilities.quadratures import gauss_legendre
from spl.linalg.stencil        import StencilVector, StencilMatrix



class InterpolationMatrix:
	"""
	
	"""

	def __init__ (self, V, coef=False, dirichlet=False):
		"""
			We are building the Interpolation Matrix
			
			
		"""
		
		if not dirichlet:
			print('InterpolationMatrix without dirichlet boundary conditions is not implemented yet')
			exit()
		# ....		
			
		
		if (coef):
			# D_i spans
			# TODO
			print('TODO')
		# ...
		else:
			# ... sizes
			[s1] = V.vector_space.starts
			[e1] = V.vector_space.ends
			[p1] = V.vector_space.pads
			
			spans_1   = V.spans
			basis_1   = V.quad_basis
			breaks = V.breaks
			
			ne = len(breaks)
			nbS = getNumberOfSplines(V)
			
			# ...
			inter = np.zeros((ne, nbS))
			# ...
			
			for ie1 in range (s1, e1 - p1 + 1):
				j_span_1 = spans_1[ie1]
				for jl_1 in range(0, p1+1):
					j = j_span_1 - p1 + jl_1
					
					inter[ie1, j] = basis_1[ie1, jl_1, 0, 0]
			
			for jl_1 in range(0, p1+1):
				inter[-1, j] = basis_1[ie1-1, jl_1, 0, 1]
			# ..
		# ...
		
#		print(inter)
#		print('\n', inter[:, 1:-1].shape)
		if (dirichlet):
			self._inter = csc_matrix(inter[:,1:-1])
		else:
			self._inter = csc_matrix(inter)
		
		self._size = inter.shape
	
	# ...
	
	def solve(self, x):
		"""
			solve the system inter*u = x
			
			Entry:
			------
				x: ndarray of shape (n, )
				
			Output:
			-------
				u : ndarray of shape (n, ) solution of the system
		"""
		
		u = spsolve(self._inter, x)
		return u
	# ..
	
	def todense(self):
		
		return self._inter.todense()
	# ...
	
	def dot (self, m):
		return self._inter.dot(m)
	# ...


# ...

class PiGrad:

	def __init__ (self, V, func, coef=False, t=0., dirichlet=False, cond=(0.,0.)):
		
		
		if not dirichlet:
			print('PiGrad without dirichlet boundary conditions is not implemented yet')
		# ...
			
		breaks = V.breaks
		
		b = np.zeros(len(breaks))
		
		# .. compute b values
		for i in range(0, len(breaks)):
			b[i] = func(breaks[i], t)
		# ...
		
		inter = InterpolationMatrix(V, coef=coef, dirichlet=dirichlet)
		# .. solve the system
		u = inter.solve(b)
		
		
		if (dirichlet):
			# we add dirichlet boundary conditions
			sol = np.zeros(u.shape[0]+2)
			
			sol[0] = cond[0]
			sol[-1] = cond[-1]
			sol[1:-1] = u
			
			self._proj = sol
		else:
			self._proj = u
		# ...	

		self._space = V
	# ...
	
	def proj (self):
		return self._proj
	# ...

	def toField(self):
		print('toField method from PiGrad class is not implemented yet')
#		return Field1d(self._proj)	
	# ....

# ...
