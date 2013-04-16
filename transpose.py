'''
Works with Matrix classes to provide low cost matrix transposition.
'''

class Transpose():

	def __init__(self,matrix):
		self.transpose = matrix
		self.T = matrix
		
	def get_data(self,i,j):
		return self.transpose._data[j][i]

	def set_data(self,i,j,val):
		self.transpose._data[j][i] = val

	def transpose(self):
		return self.T
