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


// returns the actual number of threads that Galois will use
unsigned int SPKMeansGalois::getNumThreads()
{
    return num_threads;
}


ClusterData* SPKMeansGalois::runSPKMeans(float **doc_matrix, int k, int dc, int wc)
{
    return SPKMeans::runSPKMeans(doc_matrix, k, dc, wc);
}
