from numpy import sum, sqrt, reshape, array, transpose
from matplotlib.mlab  import find
#from IPython.kernel import client

#tc = client.TaskClient()
#@tc.parallel()
def magnitude(V):
	"""	Calculate the vector magnitude of V. where V is a m x n matrix with m vectors of dimensions n.
		It assumes a matrix with row-like vectors, so that shape(V) = number of vectors, components
	"""
	#Return		
	return sqrt(dot(V,V))
	
#@tc.parallel()
def normalize(V):
	"""	Normalizes the vector V. V is a m x n matrix with m vectors of dimensions n.
		It assumes a matrix with row-like vectors, so that shape(V) = number of vectors, components
	"""
	#Return
	return V/reshape(magnitude(V),(-1,1))


#@tc.parallel()
def dot(V,B):
	"""	Dot product of V with R. V and R are m x n matrices with m vectors of dimensions n.
		It assumes matrices with row-like vectors, so that shape(V) = number of vectors, components
	"""
	#Return
	return sum(V*B,1)
	
	
def cross(A, B):
	"""Cross product of A with B. A and B are m x n matrices with m vectors of dimensions n. It assumes matrices with row-like vectors, so that shape(V) = (number of vectors, components)
	"""
	cross = array([A[:,1]*B[:,2] - A[:,2]*B[:,1], A[:,2]*B[:,0] - A[:,0]*B[:,2], A[:,0]*B[:,1] - A[:,1]*B[:,0]])
	return (transpose(cross))


#@tc.parallel()	
def GetCenterOfMass(r,coord,mass,factor):
	"""Center of mass for all cells where r > max(r)/factor
		Returns the numpy array [x,y,z]
	"""
	ind = find(r > max(r)/factor)
	x = sum(mass[ind]*coord[ind,0])/sum(mass[ind])
	y = sum(mass[ind]*coord[ind,1])/sum(mass[ind])
	z = sum(mass[ind]*coord[ind,2])/sum(mass[ind])
	return array([x,y,z])

