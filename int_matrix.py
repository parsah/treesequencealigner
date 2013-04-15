'''
A high-level class used for modeling a Matrix data-structure.
'''

from transpose import Transpose

class IntMatrix():
	''' 
	A Matrix is your traditional multi-dimensional array whereby its 
	contents are indexed by respective integers.
	'''
	def __init__(self, nrows, ncols):
		self.nrows = nrows
		self.ncols = ncols
		self.data = [[0 for _ in range(self.ncols)] 
                     for _ in range(self.nrows)]
		self.transpose = Transpose(self)
		self.T = self.transpose
    
	def get_data(self,i,j):
		return self.data[i][j]
	
	def set_data(self,i,j,val):
		self.data[i][j] = val
    
	def transpose(self):
		return self.transpose
		
