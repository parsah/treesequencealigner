'''
A high-level class used for modeling various StateMatrix data-structures.
'''

class StateMatrix():
    ''' 
    A StateMatrix is your traditional multi-dimensional array whereby its 
    contents are indexed by respective integers.
    '''
    def __init__(self, nrows, ncols):
        self.nrows = nrows
        self.ncols = ncols
        self.data = [[0.0 for _ in range(self.ncols)] 
                                for _ in range(self.nrows)]
        self.state = [[0.0 for _ in range(self.ncols)] 
                                for _ in range(self.nrows)]
        self.transpose = Transpose(self)
        self.T = self.transpose

    def get_data(self,i,j):
        ''' 
        Get the data at a specific index.
        @param i: Row number
        @param j: Column number
        '''
        return self.data[i][j]
    
    def set_data(self,i,j,val):
        self.data[i][j] = val

    def transpose(self):
        '''
        Traditional transposition of a matrix.
        '''
        return self.transpose

    def debug(self):
        ''' 
        Useful method to print-out an entire matrix row-by-row.
        '''
        for rownum in range(self.nrows):
            print(self.data[rownum], self._state[rownum])

class Transpose():
    '''
    Works with StateMatrix classes to provide low cost matrix transposition.
    '''
    
    def __init__(self,matrix):
        self.transpose = matrix
        self.T = matrix
        
    def get_data(self,i,j):
        return self.transpose.data[j][i]

    def set_data(self,i,j,val):
        self.transpose.data[j][i] = val

    def transpose(self):
        return self.T
