/* File: spkmeans_galois.cpp
 */

#include "spkmeans.h"

#include "Galois/Galois.h"
#include "Galois/Timer.h"


// constructor: set number of threads and initialize Galois
SPKMeansGalois::SPKMeansGalois(unsigned int t_) : num_threads(t_)
{
    Galois::setActiveThreads(num_threads);
    num_threads = Galois::getActiveThreads();
}



