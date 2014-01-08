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


################################################################################
# SPKMeans class
class SPKMeans():
	"""Docs"""
	
	def __init__(self, reader):
		"""Docs"""
		self.reader = reader

	def cluster(self, k):
		"""Docs"""
		self.reader.read()
		return 0
################################################################################


# define documents to be used and number of clusters
NUM_CLUSTERS = 2
docs = ["one.txt", "x", "two.txt", "three.txt"]

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