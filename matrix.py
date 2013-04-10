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
        self.data = [[0.0 for _ in range(self.ncols)] 
                     for _ in range(self.nrows)]
        self.state = [[False for _ in range(self.ncols)] 
                     for _ in range(self.nrows)]
    
    def get_data(self):
        return self.data
    
    def debug(self):
        ''' 
        Useful method to print-out an entire matrix row-by-row.
        '''
        for rownum in range(self.nrows):
            print(self.data[rownum], self.state[rownum])
            
    def transpose(self):
        ''' 
        Transpose the matrix; in other words, cell at position [2, 3] is
        now at position [3, 2]. 
        '''
        new_matrix = Matrix(nrows=self.ncols, ncols=self.nrows)
        for rownum in range(self.nrows):
            for colnum in range(self.ncols):
                # what is a row is now indexed as a column, and vice-versa.
                new_matrix.data[colnum][rownum] = self.data[rownum][colnum]
                new_matrix.state[colnum][rownum] = self.state[rownum][colnum]
        return new_matrix
        
#if __name__ == '__main__':
#    m = Matrix(nrows=5, ncols=2)
#    m.data[0][0] = 12
#    m.data[2][1] = 6
#    m.state[0][0] = True
#    m.state[2][1] = True
#    m.debug()
#    print()
#    m.transpose().debug()
    