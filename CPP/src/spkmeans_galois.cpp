/* File: spkmeans_galois.cpp
 */

#include "spkmeans.h"

#include "Galois/Galois.h"
#include "Galois/Timer.h"


// constructor: set number of threads and initialize Galois
SPKMeansGalois::SPKMeansGalois(
    float **doc_matrix_, int k_, int dc_, int wc_, unsigned int t_)
    : SPKMeans::SPKMeans(doc_matrix_, k_, dc_, wc_)
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


ClusterData* SPKMeansGalois::runSPKMeans()
{
    return SPKMeans::runSPKMeans();
}
