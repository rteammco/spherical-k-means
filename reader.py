#reader.py

import os, re


# list of filler words to be ignored
filler_words = [
	# eventually put them here
]


# global filtering regex: removes all characters matching the pattern
pattern = re.compile('[^a-z]')


class Reader():
	"""
	Reader Class:
	Used to read documents (text files) and create a set of document
	vectors, where each vector is associated with a scanned file.
	All vectors contain every seen word, including words that are not
	present in all documents. The vectors are keyed by an index that
	corresponds to a list of all seen words. The contents of the vectors
	at each index is the number of times the corresponding word has been
	seen in that document.
	
	Global variables:
	doc_vecs - list of al document vectors, filled after read() is called.
	word_list - list of all seen words, filled after read() is called.
	
	API Functions:
	__init__(files) - pass in a list of file names (paths) to be read.
	read() - read in all files in the list given to constructor.
	"""
	
	def __init__(self, files):
		"""
		Constructor:
		Parameter [files] - a list of file names (paths) to be read.
		Remembers all file names, and initializes empty lists.
		"""
		self.files = files
		self.doc_vecs = []
		self.word_list = []
	
	
	def read(self):
		"""
		Attempts to read all of the given files and compile a set of
		document vectors containing the count of every known word in
		the documents (files). A set of equal-length document vectors
		are then returned representing the occurrence count of every
		word in each document (words are identified by index only).
		Words are matched by index in the word_list list.
		"""
		for fname in self.files:
			# make sure file exists - if not, print error and skip it
			if not os.path.isfile(fname):
				print("File not found: \"{}\"".format(fname))
			else:
				with open(fname, 'r') as file:
					document_vector = []
					# read the file line by line and parse each line
					lines = file.readlines()
					for line in lines:
						words_in_line = line.split(' ')
						for word in words_in_line:
							word = self.filter(word)
							if word:
								if word in self.word_list:
									index = self.word_list.index(word)
								else:
									index = len(self.word_list)
									self.word_list.append(word)
								while len(document_vector) < index + 1:
									document_vector.append(0)
								document_vector[index] += 1
								#print(document_vector)
					self.doc_vecs.append(document_vector)
		
		# pad the vectors so that they all have equal length
		self.pad()
	
	
	def filter(self, word):
		"""
		Converts all words into lowercase form, and removes all excluded
		characters as defined by the global "pattern" regular expression.
		If a word is an empty string or if matched as a filler word (as
		defined by the "filler_words" list), "None" will be returned.
		Otherwise, the modified word will be returned.
		"""
		word = word.strip() # remove all whitespace
		word = word.lower() # make word all lowercase
		word = pattern.sub('', word) # keep only lowercase letters
		if not word or word in filler_words:
			word = None
		return word
	
	
	def pad(self):
		"""
		Pad the space after each document vector such that all document
		vectors end up having the exact same number of words. Those
		words in the end will be set to count 0.
		"""
		for document_vector in self.doc_vecs:
			while len(document_vector) < len(self.word_list):
				document_vector.append(0)