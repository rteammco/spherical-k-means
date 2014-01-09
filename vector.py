import math


class Vector():
	
	
	def __init__(self, size):
		self.size = size
		self.array = [0] * size
	
	
	def __len__(self):
		return self.size
	
	
	def __add__(self, other):
		result = list(self.array)
		for i in range(self.size):
			result[i] += other[i]
		return result
	
	
	def __mul__(self, scalar):
		result = list(self.array)
		for i in range(self.size):
			result[i] *= scalar
		return result
	
	
	def norm(self):
		result = list(self.array)
		n = 0
		for item in self.array:
			n += item * item
		n = math.sqrt(n)
		for i in range(self.size):
			result[i] /= n
		return result