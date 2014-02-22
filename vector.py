# Vector class:
# Supports basic mathematical vector operations as needed by the
# implementing algorithms. This is NOT a complete vector math library.

import math


class Vector():
	"""
	Vector class supports basic vector functionality, defined by the
	following functions and operators; "v" is the Vector:
	
	append: 	insert a new element to the end of this Vector.
	v[i]:		set or get item at index i in this Vector.
	len(v):		returns the length (size) of this Vector.
	v + v2:		adds this Vector to Vector v2 and returns the result.
	v * s:		multiplies this Vector by scalar s.
	dot(b):		returns a dot product of this Vector with Vector b.
	norm:		returns a normalized version of this Vector.
	print(v):	prints a string representation of this Vector.
	"""
	
	
	def __init__(self, size=0, val=0):
		"""
		Initialize a zero Vector of the given size. Vector can be
		initialized with a non-zero value if "val" parameter is set.
		If size is not provided, vector will initially be empty.
		"""
		self.size = size
		if size > 0:
			self.array = [val] * size
		else:
			self.array = []
		
		# optimizations: if marked as modified (true), norm and other values
		#	must be re-calculated when used.
		self.modified = True
		self.norm_val = 0
	
	
	@staticmethod
	def from_array(array):
		"""STATIC: construct a new Vector from the given array."""
		new_vec = Vector()
		new_vec.set_array(array)
		return new_vec
	
	
	def set_array(self, array):
		"""Set the contents of this Vector to the given array."""
		self.array = array
		self.size = len(self.array)
		self.modified = True
	
	
	def append(self, val):
		"""Inserts "val" to the end of this vector. Increments size."""
		self.array.append(val)
		self.size += 1
		self.modified = True
	
	
	def __getitem__(self, index): # [] get operator
		"""
		Returns the item in this Vector at the given index.
		Causes an regular list error if index is out of bounds.
		"""
		return self.array[index]
	
	
	def __setitem__(self, index, val): # [] set operator
		"""
		Sets the item in this Vector at the given index to the given value.
		Causes a regular list error if index is out of bounds.
		"""
		self.array[index] = val
		self.modified = True
		
	
	def __len__(self): # len() operator
		"""Returns length (size) of this Vector."""
		return self.size
	
	
	def __add__(self, other): # + operator
		"""
		Adds Vector "other" to this Vector and returns the result.
		Both vectors are expected to be the same size (length), otherwise
		an AssertionError will be thrown.
		"""
		assert other.size == self.size
		result = Vector.from_array(self.array)
		for i in range(self.size):
			result[i] += other[i]
		self.modified = True
		return result
	__radd__ = __add__
	
	
	def __mul__(self, scalar): # * operator
		"""Multiplies this Vector with the scalar and returns the result."""
		result = Vector.from_array(self.array)
		for i in range(self.size):
			result[i] *= scalar
		return result
	__rmul__ = __mul__
	
	
	def dot(self, other):
		"""
		Returns the dot product NUMBER value of this Vector with "other".
		Both vectors are expected to be the same size (length), otherwise
		an AssertionError will be thrown.
		"""
		assert other.size == self.size
		result = 0
		for i in range(self.size):
			result += self.array[i] * other[i]
		return result
		
	
	def norm(self):
		"""Returns the norm of this Vector."""
		#result = Vector.from_array(self.array)
		n = 0
		for item in self.array:
			n += item * item
		n = math.sqrt(n)
		#for i in range(self.size):
		#	result[i] /= n
		#return result
		return n
	
	
	def normalize(self):
		"""
		Sets this vector to a normalized version of itself.
		WARNING: All old data is lost after running this function.
		"""
		n = self.norm()
		for i in range(self.size):
			self.array[i] /= n
		
	
	def __str__(self):
		"""Returns a string representation of this Vector."""
		return self.array.__str__()
