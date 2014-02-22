#!/usr/bin/env python3

################################################################################
# A Python implementation of the Spherical K-Means clustering algorithm.	
#
# The SPK-Means algorithm is used to create K clusters into which it sorts
# given documents containing text data. These clusters are generated as a
# means of finding patterns in large amounts of data.
# The number of clusters, K, is provided as a parameter.
#
# Implementation by Richard Teammco
# teammco@cs.utexas.edu
#
# Tested with Python versions: 3.1.3, 3.2.3
#
# Doctest:
#   $ python3 -m doctest spkmeans.py
################################################################################


# import Reader module
from reader import Reader
from vector import Vector

# import global modules
import math


################################################################################
# SPKMeans class
class SPKMeans():
	"""
	Contains functions to cluster a documenting using the Spherical
	K-Means algorithm.
	"""
	
	
	def __init__(self, reader):
		"""Sets up the reader and initializes class variables."""
		self.reader = reader
		self.doc_vecs = None
		self.num_docs = 0
		self.num_words = 0
		

	def read_docs(self):
		"""Reads the documents from files into memory, and report stats."""
		self.doc_vecs = self.reader.read()
		self.num_docs = len(self.doc_vecs)
		self.num_words = len(self.reader.word_list)
		self.reader.report()
		
		
	def cluster(self, k):
		"""Run the full SPKMeans algorithm, and return the partitions."""
		
		# k must be at least 2 (otherwise meaningless)
		if k < 2:
			print("Warning: must use at least 2 partitions. Stopping.")
			return None
		
		self.read_docs()
		
		# number of documents must be at least 2 (otherwise meaningless)
		if self.num_docs < 2:
			print("Warning: must use at least 2 documents. Stopping")
			return None
		
		print("Running SPKMeans clustering: {} partitions.".format(k))
		
		# normalize the document vectors
		for doc_v in self.doc_vecs:
			doc_v.normalize()
		
		# TODO - txn scheme
		
		# create first partition set and concept vectors
		partitions = self.randomize_partitions(k)
		concepts = self.get_concepts(partitions)
		
		while True: # while delta_quality < q_threshold:
			partitions = [[]] * k
			
			for i in range(self.num_docs): # for each document
				#find closest doc vector of cosine similarity
				closest = -1
				closest_val = 0
				for cv_index in range(k):
					scrap = concepts[cv_index]
					dot_p = self.doc_vecs[i].dot(concepts[cv_index])
					if closest == -1 or dot_p < closest_val:
						closest = cv_index
						closest_val = dot_p
						
				# put the document Vector into the partition associated
				#	with the closest concept Vector
				partitions[closest].append(self.doc_vecs[i])
				
			# recalculate concept vectors
			concepts = self.get_concepts(partitions)
			return 123 # TODO - remove
			
			#
			# NOTE: If a partition ends up being empty, there can be no
			#	concept vector for it (norm of 0 vector is undefined)...
			#	Can an empty partition happen? If so, how to deal with it?
			#
		
		print(concepts)
		return 0
	
	
	def txn(self, doc_vec_i):
		"""
		Scales the given document Vector using the TXN scheme.
		doc_vec_i[j] = t_ji * g_j * s_i
			- (doc i, word j)
		t_ji:	term weighting component
			- depends on count of word j in doc i
		g_j:	global weighting component
			- depends on number of docs containing word j
		s_i:	normalization component for doc i
			= sqrt[sum(j,word_count):(t_ji * g_j)^2]
		
		*** TESTING: ***
		###### Set up class:
		>>> s = SPKMeans(None)
		>>> s.num_words = 5
		
		###### Test vectors:
		>>> v = [1,2,3,0,0]
		>>> s.txn(v)
		>>> print(v)
		[0,0,0,0]
		>>> v = [0,0,2,3,4]
		>>> s.txn(v)
		>>> print(v)
		[0,0,0,0]
		"""
		# calculate s_i (normalization component) for vector (i):
		#	= sqrt[sum(j,word_count):(t_ji * g_j)^2]
		sum_v = 0
		for j in range(self.num_words):
			t_ji = doc_vec_i[j] # term weighting component
			g_j = 1 # global weighting component (TODO - ignored)
			sum_v += (t_ji*g_j)*(t_ji*g_j) # = math.pow((t_ji*g_j), 2)
		s_i = math.sqrt(sum_v)
		print(s_i)
		
		# calculate the new vector values for vector (i):
		#	doc_vec_i[j] = t_ji * g_j * s_i
		for j in range(self.num_words):
			t_ji = doc_vec_i[j]
			g_j = 1
			doc_vec_i[j] = t_ji * g_j * s_i
		
		
	def get_concepts(self, partitions):
		"""Returns a set of concept Vectors, one for each partition."""
		concepts = []
		for p in partitions:
			concepts.append(self.compute_concept(p))
		return concepts
	
	
	def compute_concept(self, p):
		"""
		Computes the Concept Vector of given partition p, and returns said
		Concept Vector. Parameter p must be a list containing at least one
		document Vector.
		"""
		# computer sum of all vectors in partition p
		cv = Vector(self.num_words)
		for doc_v in p:
			cv += doc_v
			
		# compute the mean vector for partition using the sum vector
		cv *= (1/len(p))
		
		# computer the norm of the mean vector
		cv.normalize()
		
		return cv
	
	
	def randomize_partitions(self, k):
		"""
		Create a set of k partitions with the documents randomly
		distributed to a cluster. Maintains as much distribution equality
		between the partitions as possible.
		Returns the set of partitions.
		TODO: not random; splits up by order docs were scanned in
		"""
		partitions = []
		num_per_partition = int(self.num_docs / k)
		for pindex in range(k):
			partition_vecs = []
			partition_size = num_per_partition
			if pindex == k-1:
				partition_size = partition_size + (self.num_docs % k)
			for num in range(partition_size):
				num += (pindex * (k-1))
				partition_vecs.append(self.doc_vecs[num])
			partitions.append(partition_vecs)
		
		return partitions
################################################################################


# start here
def main():
	# define documents to be used and number of clusters
	# TODO - docs and number of clusters as parameters
	NUM_CLUSTERS = 2
	docs = ['one.txt', 'x', 'two.txt', 'three.txt']
			#'long1.txt', 'long2.txt', 'long3.txt']
	
	# set up reader and run SPKMeans clustering
	reader = Reader(docs)
	clusters = SPKMeans(reader).cluster(NUM_CLUSTERS)
	print(clusters)


main()


#
