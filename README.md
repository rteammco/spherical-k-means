Spherical K-Means Algorithm
=======

An implementation of the spherical K-means algorithm in Matlab and C++.

The C++ version emphasizes a multithreaded implementation, and features three ways of running the algorithm. It can be executed with a single-thread (same as the Matlab implementation), or using OpenMP or Galois (http://iss.ices.utexas.edu/?p=projects/galois). The purpose of this code is to optimize and compare the different parallel paradigms to maximize the efficiency of the algorithm.

See the subdirectories for specific instructions on using those programs and data sets.

The basic version of K-means was implemented as described in this paper: http://www.cs.utexas.edu/users/inderjit/public_papers/concept_mlj.pdf
