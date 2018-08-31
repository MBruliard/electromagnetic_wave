
import numpy as np


def listToVectorColumn(liste):

	"""
		Change a list into a vector Column 
	"""
	
	size = len(liste)
	vector = np.asarray(liste)
	return vector
# --

# -- 
def vectorToList(vectorColumn):
	"""
		Change a vector into a list
	"""
#	vectorLine = vectorColumn.transpose()
	liste = vectorColumn.tolist()
	
	return liste
# --

# -- 
def addDirichletConditions(vector, leftdirichlet, rightdirichlet):
	"""
		Requires:
		--------
			* vector is a np.array 2D of size (n, 1)
			* leftdirichlet is a real
			* rightdirichlet is a real
		
		Returns:
		--------
			* newVector is a np.array 2D of size (n+2, 1) where newVector(1,1) = leftdirichlet and newVector(n+2, 1) = rightdirichlet
	"""
	
	newvector = vectorToList(vector)
	newvector.insert(0, leftdirichlet)
	newvector.append(rightdirichlet)
	newvector = listToVectorColumn(newvector)
	
	return newvector
# --

# -- 
def supprDirichletConditions(vector):
	"""
		Requires:
		--------
			* vector is a np.array 2D of size (n, 1)
		
		Returns:
		--------
			* newVector is a np.array 2D of size (n-2, 1) 
	"""
	
	return vector[1:-1]
# --


# --
def group_two_vectors(v, w):
	"""
		Requires:
		---------
			v is a np.ndarray of size (n, 1)
			w is a np.ndarray of size (m, 1)
		
		Returns:
		--------
			res is a np.ndarray of size (n+m, 1)
	"""
	
	n = v.shape[0]
	m = w.shape[0]
	
	res = np.zeros(n+m)
	
	res[0:n] = v
	res[n:m+n] = w
	
	return res;
# --


# --	
def ungroup_vector(v, separator_index):
	"""
		Requires:
		---------
			v is a np.ndarray 2D of size (n, 1)
			separator_index is int between {1, ..., n-1}	
		
		Returns:
		--------
			fisrt is a np.ndarray 2D of size (separator_index, 1)
			second is a np.ndarray 2D of size (n-separator_index, 1)
	"""
	n = v.shape[0]
	first = np.zeros(separator_index)
	second = np.zeros(n-separator_index)
	
	first = v[0:separator_index]
	second = v[separator_index:n]
	
	
	return first, second
# --















