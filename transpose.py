'''
Works with Matrix classes to provide low cost matrix transposition.
'''

class Transpose():

	def __init__(self,matrix):
		self.transpose = matrix
		self.T = matrix
		
	def getData(self,i,j):
		return self.transpose.data[j][i]

	def setData(self,i,j,val):
		self.transpose.data[j][i] = val

	def transpose(self):
		return self.T
