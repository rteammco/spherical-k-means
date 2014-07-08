/* File: spkmeans_openmp.cpp
 *
 * Contains the definitions of the OpenMP version of the spherical K-means
 * algorithm. OpenMP provides multithreading.
 */

#include "spkmeans.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include <omp.h>
#include "Galois/Timer.h"

#include "cluster_data.h"
#include "vectors.h"

using namespace std;


// CONSTRUCTOR: set a pre-defined number of threads.
SPKMeansOpenMP::SPKMeansOpenMP(
    float **doc_matrix_, int k_, int dc_, int wc_, unsigned int t_)
    : SPKMeans::SPKMeans(doc_matrix_, k_, dc_, wc_)
{
    // make sure num_threads doesn't exceed the max available
    if(t_ > omp_get_max_threads())
        num_threads = omp_get_max_threads();
    // and it can't be less than 1 (if <= 0, set to max)
    else if(t_ <= 0)
        num_threads = omp_get_max_threads();
    else
        num_threads = t_;
    omp_set_num_threads(num_threads);
}


// returns the actual number of threads that OpenMP will use
unsigned int SPKMeansOpenMP::getNumThreads()
{
    return num_threads;
}


// Overridden for parallelizing:
/*float SPKMeansOpenMP::computeQ(
    float ***partitions, int *p_sizes, float **concepts)
{
    float quality = 0;
    //mutex mut;
    //#pragma omp parallel for num_threads(num_threads)
    for(int i=0; i<k; i++) {
      //  mut.lock();
        quality +=
            SPKMeans::computeQ(partitions[i], p_sizes[i], concepts[i]);
      //  mut.unlock();
    }
    return quality;
}*/


// Runs the spherical k-means algorithm on the given sparse matrix D and
// clusters the data into k partitions.
ClusterData* SPKMeansOpenMP::runSPKMeans()
{
    // keep track of the run time for this algorithm
    Galois::Timer timer;
    timer.start();

    // keep track of all individual component times for analysis
    Galois::Timer ptimer;
    Galois::Timer ctimer;
    Galois::Timer qtimer;
    float p_time = 0;
    float c_time = 0;
    float q_time = 0;

    // apply the TXN scheme on the document vectors (normalize them)
    txnScheme();

    // initialize the data arrays; keep track of the arrays locally
    ClusterData *data = new ClusterData(k, dc, wc, doc_matrix);
    float **concepts = data->concepts;
    bool *changed = data->changed;
    float *cValues = data->cValues;
    float *qualities = data->qualities;

    // compute initial partitions, concepts, and quality
    initPartitions(data);
    float quality = computeQ(data);
    cout << "Initial quality: " << quality << endl;


    // do spherical k-means loop
    float dQ = Q_THRESHOLD * 10;
    int iterations = 0;
    while(dQ > Q_THRESHOLD) {
        iterations++;

        // compute new partitions based on old concept vectors
        ptimer.start();
        #pragma omp parallel for
        for(int i=0; i<dc; i++) {
            int pIndx = 0;
            // only update cosine similarities if partitions have changed
            if(changed[0])
                cValues[i*k] = cosineSimilarity(concepts[0], i);
            for(int j=1; j<k; j++) {
                if(changed[j]) // again, only if changed
                    cValues[i*k + j] = cosineSimilarity(concepts[j], i);
                if(cValues[i*k + j] > cValues[i*k + pIndx])
                    pIndx = j;
            }
            data->assignPartition(i, pIndx);
        }
        ptimer.stop();
        p_time += ptimer.get();

        // update which partitions changed since last time, then swap pointers
        data->findChangedPartitions();
        data->swapAssignments();

        // compute new concept vectors
        ctimer.start();
        #pragma omp parallel for
        for(int i=0; i<k; i++) {
            // only update concept vectors if partition has changed
            if(changed[i]) {
                delete[] concepts[i];
                concepts[i] = computeConcept(data, i);
            }
        }
        ctimer.stop();
        c_time += ctimer.get();

        // compute quality of new partitioning
        qtimer.start();
        float n_quality = computeQ(data);
        dQ = n_quality - quality;
        quality = n_quality;
        qtimer.stop();
        q_time += qtimer.get();

        // report the quality of the current partition
        reportQuality(data, quality, dQ);
    }


    // report runtime statistics
    timer.stop();
    reportTime(iterations, timer.get(), p_time, c_time, q_time);

    // return the resulting partitions and concepts in the ClusterData struct
    return data;
}
