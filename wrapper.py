from matrix import Matrix

class DirectionalMatrixWrapper():
	'''
	A wrapper for directional score matrices
	'''
	
	def __init__(self,nrows,ncols,T=None):
		self.nrows = nrows
		self.ncols = ncols
		if T is None:
			self.score = Matrix(nrows,ncols)
			self.extend_flag = Matrix(nrows,ncols)
			self.T = DirectionalMatrixWrapper(nrows,ncols,self)
		else:
			self.T = T
			self.score = T.score.T
			self.extend_flag = T.extend_flag.T

	def transpose(self):
		return self.T
