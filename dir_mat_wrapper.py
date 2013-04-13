'''
A wrapper for directional score matrices
'''

class DirMatWrapper():

	def __init__(self,nrows,ncols,T=None):
		self.nrows = nrows
		self.ncols = ncols
		if T is None:
			self.score = IntMatrix(nrows,ncols)
			self.flag = IntMatrix(nrows,ncols)
			self.T = DirMatWrapper(nrows,ncols,self)
		else:
			self.T = T
			self.score = T.score.T
			self.flag = T.flag.T

	def transpose(self):
		return self.T
		
