Spherical K-Means: C++ Version
=======

This is a C++ implementation of the spherical K-means algorithm.


Dependencies
-------

**OpenMP** (*required*) - available by using a newer version of gcc/g++ to compile the code. Most Linux distributions will support this library.

**Boost** (*required*) - needed by the `Timer` object (`src/timer.h/cpp`). If you don't want to use boost, you can change the Timer class as long as it conforms to the public specifications defined in `timer.h`.

**Galois** (*optional*) - if you do not have this library installed, the Makefile will build the code without it (see below), but you will not be able to run the Galois version. Galois is an open source project available for download here: http://iss.ices.utexas.edu/?p=projects/galois. There is no performance change when using Galois as opposed to OpenMP in this implementation (hence why I started the Galois-only version). NOTE: Installing Galois requires a newer version of cmake.


Build
-------

To build the code, you generally only need to run `make`. If you choose to install Galois, set `GALOIS_PATH_RAW` to the appropriate directory on line 5 of the Makefile. NOTE: On a 64-bit system, the Makefile will append "_64" to that path name. You might want to simply delete lines 5-13, and set `GALOIS_PATH` to wherever you installed Galois. If Galois is not installed, the Makefile will automatically ignore it and build the code without it.

The compiler is specified in the Makefile on line 18, and all required flags are specified on line 19.


Running the Code
-------

All runtime flags can be viewed by running `./spkmeans --help`. In general, you will probably want to run it with the following options:

`./spkmeans -d path/to/docfile -k n --noresults`

Here, `path/to/docfile` is the input document data. For example, using the provided data sets, you can use `../TestData/documents`.

`n` is the size of k (i.e. number of clusters). You can instead use the `--autok` flag which will try to approximate the optimal number of clusters given the data set.

If `--noresults` is not provided, the program will print out the top 10 words for each resulting cluster. If this option *is* set, the program will still print general clustering statistics.

Adding a vocabulary file with `-v path/to/vocabfile` is just useful if don't silence the results. Instead of just printing word IDs, it will print the actual words themselves. For the provided data sets, only one vocabulary file is given (`TestData/vocabulary`), associated with the `TestData/documents` data set.

If you want to **run with multiple threads**, add either the `--openmp` or `--galois` flags. If these aren't provided, the program will automatically run the single-threaded version. By default, both methods will use the maximum number of threads available. To specify a different number of threads, add runtime option `-t n` where `n` is the number of threads.

All other options are fairly unimportant. `--noscheme` will skip the normalization step for data sets that are already normalized. If a data set has already been normalized, it will just waste a few seconds of time. `--noop` will disable all optimizations in the algorithm. This will obviously make it slower, and does not impact the results. The flag exists for testing purposes only.


Source Code
-------

All original source code is in the `src` directory. The README file in that directory provides basic documentation on the organization of the source files.
