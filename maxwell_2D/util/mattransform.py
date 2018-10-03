
import sys
import numpy as np

def fourToTwoD (mat):
	"""
		mat est un ndarray 4D
	"""
	(m, n, p, q) = mat.shape
	
	mmax = m*p
	nmax = n*q
	
	res = np.zeros((mmax, nmax))

	for i1 in range(0, m):
		for i2 in range(0, n):
			i = i1*p
			j = i2*q
			res[i:i+p, j:j+q] = mat[i1, i2, :, :]
	# ...
				
	return res
# ...


	
def twoToOneD (mat):
	"""
		Transform a 2D matrix into a vector column:
		
		Entry:
		------
			mat: 2D ndarray of shape (m, n)
			
		Ouput:
		-------
			vect: 1D ndarray of shape (m*n, )
	"""


	(m, n) = mat.shape
	
	vect = np.zeros(m*n)
	
	for i in range(0, m):
		for j in range(0, n):
			vect[i*n+j] = mat[i, j]
	# ...
	
	return vect
# ....	

	
def nonameyet (mat, field):
	"""
		
		
		Entries:
		--------
			mat: ndarray 4D of shape (m, n, p, q)
			
			field: ndarray of shape (p, q)
			
		Output:
		-------
			y: ndarray of shape (m, n)
				y[i, j] = sum_p sum_q mat[i, j, p, q]*field[p, q]
			
	"""
	(m, n , p, q) = mat.shape
	(fp, fq) = field.shape
	
	if (not (fp==p and fq==q)):
		print("ERROR in nonameyet: dimension not match")
		sys.exit()
	
	y = np.zeros((m, n))
	
	for i1 in range(0, m):
		for i2 in range(0, n):
			
			s = 0.
			for j1 in range(0, p):
				for j2 in range(0, q):
					s = s + mat[i1, i2, j1, j2]*field[j1, j2]
			y[i1, i2] = s 
	# ...
	return y
	
	
	
