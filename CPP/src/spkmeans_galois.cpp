/* File: spkmeans_galois.cpp
 */

#include "spkmeans.h"

#include "Galois/Galois.h"
#include "Galois/Timer.h"


// constructor: set number of threads and initialize Galois
SPKMeansGalois::SPKMeansGalois(unsigned int t_)
{
    // if number of threads given is <= 0, set to max
    if(t_ <= 0)
        t_ = 2; // TODO - no support in API for this?
    Galois::setActiveThreads(t_);
    num_threads = Galois::getActiveThreads();
}


// returns the actual number of threads that Galois will use
unsigned int SPKMeansGalois::getNumThreads()
{
    return num_threads;
}


ClusterData* SPKMeansGalois::runSPKMeans(float **doc_matrix, int k, int dc, int wc)
{
    return SPKMeans::runSPKMeans(doc_matrix, k, dc, wc);
}
