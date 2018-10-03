__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-25"
__aim__ = "Field class for maxwell 2D"


import numpy as np

from spl.fem.splines import SplineSpace
from spl.linalg.stencil import StencilVector, StencilMatrix
from spl.fem.tensor  import TensorFemSpace


class Field1d:

#	def __init__(self, V, func, t):
#		
#		self._values = valuesfromgreville(V, func, t=t)
#	# ...
	
	def __init__ (self, values):
#		"""
#			Entries:
#			--------
#				values: ndarray of shape (n,)
#		"""
#	
#		self._values = values

		"""
		
			Entries: 
			--------
				values: ndarray of shape (n , m)
						E[j1, j2] = e_{j1, j2}
		"""
		
		self._values = values
		

	# ...
	
	def dir(self):
		return self._values
	# ...
	
	def size(self):
		return self._values.shape
	# ...
	
#	def tovectorColumn(self):
#		(m, n) = self.size()
#		val = np.zeros(m*n)
#		
#		for i in range(0, m):
#			for j in range(0, n):
#				val[i*n+j] = self._values[i, j]
#				
#		return val
#	# ...

# ...


class Field2d:

#	def __init__(self, V1, V2 , func_1, func_2, t=0.):
#	
#		val1 = valuesfromgreville(V1, func_1, t=t)
#		val2 = valuesfromgreville(V2, func_2, t=t)
#		
#		self._values = [val1, val2]
#	# ...
	
	
	def __init__ (self, val1, val2):

		"""
			Entries:
			--------
				val1: np.ndarray of shape (n ,)
				
				val2 : np.ndarray of shape (p, )
		"""
		self._values = [val1, val2]
	# ...
	
	def dirx (self):
		return self._values[0]
	# ...
	
	def diry (self):
		return self._values[-1]
	# ...
	
	def size(self):
		return (self._values[0].shape, self._values[-1].shape)
	
#	def tovectorColumn(self):
#		((m, n), (p, q)) = self.size()
#		
#		_val1 = self._values[0]
#		_val2 = self._values[1]
#		
#		val1 = np.zeros(m*n)
#		val2 = np.zeros(p*q)
#		
#		for i in range(0, m):
#			for j in range(0, n):
#				val1[i*n+j] = _val1[i, j]
#		
#		for i in range(0, p):
#			for j in range(0, q):
#				val2[i*q+j] = _val2[i, j]
#			
#			
#		return [val1, val2]
#	# ...
# ... 


def valuesfromgreville (V, func, t=0. ):
	"""
		Compute y_i = func(x_i, t) for x_i greville abscissae
		
		Entries:
		---------
			V: SplineSpace or TensorFemSpace
			
			func: function studing
			
			t: value of the time variable
			
		Output:
		-------
			y: ndarray of shape(n, ) 
	"""
	
	[V1, v2] = [W for W in V.spaces]
	[greville_1, greville_2] = [W.greville for W in V.spaces]
	n1 = len(greville_1)
	n2 = len(greville_2)
	
#	values = np.ndarray(n1*n2)
#	
#	for i in range(0, n1):
#		for j in range(0, n2):
#			xi = greville_1[i]
#			yj = greville_2[j]
#			
#			values[i*n2 + j] = func(xi, yj, t)
#	# ...
	
	values = np.zeros((n1, n2))
	
	for i1 in range (n1):
		for i2 in range(n2):
			x = greville_1[i1]
			y = greville_2[i2]
			
			values[i1, i2] = func(x, y, t)
	# ...
	
	return values
# ...

