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
# Tested with Python version 3.1.3
################################################################################


# import Reader module
from reader import Reader
from vector import Vector


################################################################################
# SPKMeans class
class SPKMeans():
	"""Docs"""
	
	
	def __init__(self, reader):
		"""Docs"""
		self.reader = reader
		self.reader.read()
		self.reader.report()
		
		self.partitions = []
		

	def cluster(self, k):
		"""Docs"""
		self.randomize_partitions(k)
		return 0
	
	
	def randomize_partitions(self, k):
		"""
		Create a set of k partitions with the documents randomly
		distributed to a cluster. Maintains as much distribution equality
		between the partitions as possible.
		TODO: not random; splits up by order docs were scanned in
		"""
		num_per_partition = int(self.reader.num_docs / k)
		for pindex in range(k):
			partition_vecs = []
			partition_size = num_per_partition
			if pindex == k-1:
				partition_size = partition_size + (self.reader.num_docs % k)
			for num in range(partition_size):
				num += (pindex * (k-1))
				partition_vecs.append(self.reader.doc_vecs[num])
			self.partitions.append(partition_vecs)
		#for p in self.partitions:
		#	print(p)
	
	
	def compute_concept(self, p):
		"""
		Computers the concept vector of given partition p, and returns said
		concept vector. Parameter p must be a list containing at least one
		document vector (list).
		"""
		p_size = len(p)
		num_words = len(self.reader.word_list)
		
		# computer sum of all vectors in partition p
		sum_v = [0] * num_words # list starts with all zeros
		for doc_v in p:
			for w in range(num_words):
				sum_v[w] += doc_v[w]
		
		# compute the mean vector for partition using the sum vector
		mean_v = [1 / p_size] * num_words
		for w in range(num_words):
			mean_v[w] *= sum_v[w]
		
		# computer the norm of the mean vector
		norm_v = []
		cv = 0
		return cv
################################################################################


# define documents to be used and number of clusters
NUM_CLUSTERS = 2
docs = ['one.txt', 'x', 'two.txt', 'three.txt']
		#'long1.txt', 'long2.txt', 'long3.txt']

# start here
if __name__ == "__main__":
	reader = Reader(docs)
	clusters = SPKMeans(reader).cluster(NUM_CLUSTERS)
	print(clusters)
	

#	(1) Start with arbitrary partitioning 0: {ℼ(0)j}kj=1.
#		{c(0)j}kj=1 denotes concept vectors associated with given partitioning 0.
#		Set t = 0 (index of iteration).
#	(2) ∀xi (1 to n), find concept vector closest in cosine similarity to xi.
#		Compute new partitioning {ℼ(t+1)j}kj=1 induced by old concept vectors {c(t)j}kj=1.
#		“ℼ(t+1)j is the set of all document vectors closest to the concept vector c(t)j.”
#		- if doc vector closest to more than one, assign it randomly to one of the clusters.
#	(3) Compute new concept vectors for (t+1):
#		c(t+1)j = m(t+1)j / || m(t+1)j ||, 1 ≤ j ≤ k.						( 8 )
#		m(t+1)j = centroid or mean of document vectors in cluster ℼ(t+1)j.
#	(4) 	If stopping criterion met: set ℼ*j = ℼ(t+1)j and c*j = c(t+1)j (for all 1 ≤ j ≤ k). EXIT.
#		Else: t++; GOTO step_2.
#		stopping criterion may be difference in Q at t and t+1 <= threshold.