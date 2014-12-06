Spherical K-Means Algorithm
=======

An implementation of the spherical K-means algorithm in Matlab and C++.

The C++ version emphasizes a multithreaded implementation, and features three ways of running the algorithm. It can be executed with a single-thread (same as the Matlab implementation), or using OpenMP or Galois (http://iss.ices.utexas.edu/?p=projects/galois). The purpose of this code is to optimize and compare the different parallel paradigms to maximize the efficiency of the algorithm.

See the subdirectories for specific instructions on using those programs and data sets.

The basic version of K-means was implemented as described in this paper: http://www.cs.utexas.edu/users/inderjit/public_papers/concept_mlj.pdf

Directories:
-------

`CPP` contains the C++ version of the code as well as documentation on installing and running it. Galois is optional for running this code.

`Galois` contains an experimental Galois-only version of the algorithm. This implementation is not finished, and requires Galois to work.

`Matlab` contains a prototype Matlab version of the algorithm. This version is not optimized, and has fairly poor performance.

`Notes` can be ignored. It is just a folder where I put my personal notes when developing this code.

`Scripts` contains scripts to automate experiments and performance tests of the code.

`TestData` contains various input data sets for testing. The `*.gr` files are required for the Galois version. All other files have been formatted to work with the CPP version.
