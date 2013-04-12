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
    
		return self.transpose
	def transpose(self):
	
    def debug(self):
        ''' 
        Useful method to print-out an entire matrix row-by-row.
        '''
        for rownum in range(self.nrows):
            print(self.data[rownum], self.state[rownum])
            
	
