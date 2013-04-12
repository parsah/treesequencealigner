'''
A wrapper for directional score matrices
'''

class DirMatWrapper():

  def __init__(self,nrows,ncols):
		self.nrows = nrows
		self.ncols = ncols
		self.score = IntMatrix(nrows,ncols)
		self.flag = IntMatrix(nrows,ncols)
