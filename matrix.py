'''
A high-level class used for modeling various StateMatrix data-structures.
@author: Todd Gillette and Parsa Hosseini
'''

class AbstractMatrix():
    """
    A high-level matrix with capability of setting row and column names and
    dimensional information.
    """

    def __init__(self, nrows, ncols):
        """
        An abstract matrix capable of having row and column data points.
        """
        self.nrows = nrows
        self.ncols = ncols
        self.data = [[0.0 for _ in range(self.ncols)] 
                     for _ in range(self.nrows)]

    def get_width(self):
        """
        Retrieve how many columns there are in the current matrix.
        @return: matrix width (number of columns).
        """
        if len(self.data) == 0:
            return 0
        else:
            return len(self.data[0])

    def get_height(self):
        """
        Retrieve the number of rows comprising the current matrix.
        @return: matrix height (number of rows).
        """
        return len(self.data)

    def is_square(self):
        """
        Determines if all rows in the matrix are the same length.
        @return: boolean.
        """
        return len(set([len(row) for row in self.data])) == 1

class StateMatrix(AbstractMatrix):
    ''' 
    A StateMatrix is your traditional multi-dimensional array whereby its 
    contents are indexed by respective integers.
    '''
    def __init__(self, nrows, ncols):
        super(StateMatrix, self).__init__(nrows, ncols)
        self.state = [[0.0 for _ in range(self.ncols)] 
                                for _ in range(self.nrows)]
        self.transpose = TranspositionFactory(self)
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

class TranspositionFactory():
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

class DirectionalMatrixWrapper():
    '''
    A wrapper for directional score matrices
    '''
    
    def __init__(self,nrows,ncols,T=None):
        self.nrows = nrows
        self.ncols = ncols
        if T is None:
            self.score = StateMatrix(nrows,ncols)
            self.extend_flag = StateMatrix(nrows,ncols)
            self.T = DirectionalMatrixWrapper(nrows,ncols,self)
        else:
            self.T = T
            self.score = T.score.T
            self.extend_flag = T.extend_flag.T

    def transpose(self):
        return self.T

