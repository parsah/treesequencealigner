'''
A high-level class used for modeling a Matrix data-structure.
'''

class Matrix():
	''' 
	A Matrix is your traditional multi-dimensional array whereby its 
	contents are indexed by respective integers.
	'''
	def __init__(self, nrows, ncols):
		self.nrows = nrows
		self.ncols = ncols
		self._data = [[0.0 for _ in range(self.ncols)] 
                     for _ in range(self.nrows)]
    
	def get_data(self,i,j):
		return self._data[i][j]
	
	def set_data(self,i,j,val):
		self._data[i][j] = val
    
	def transpose(self):
		return self.transpose
    
	def debug(self):
		''' 
		Useful method to print-out an entire matrix row-by-row.
		'''
		for rownum in range(self.nrows):
			print(self._data[rownum], self._state[rownum])

#if __name__ == '__main__':
#    m = Matrix(nrows=5, ncols=2)
#    m.data[0][0] = 12
#    m.data[2][1] = 6
#    m.state[0][0] = True
#    m.state[2][1] = True
#    m.debug()
#    print()
#    m.transpose().debug()
    
