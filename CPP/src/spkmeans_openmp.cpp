/* File: spkmeans_openmp.cpp
 *
 * Contains the definitions of the OpenMP version of the spherical K-means
 * algorithm. OpenMP provides multithreading.
 */

#include "spkmeans.h"

#include <iostream>

#include <omp.h>
#include "timer.h"

#include "cluster_data.h"

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



// Runs the spherical k-means algorithm on the given sparse matrix D and
// clusters the data into k clusters.
ClusterData* SPKMeansOpenMP::runSPKMeans()
{
    // keep track of the run time for this algorithm
    Timer timer;
    timer.start();

    // keep track of all individual component times for analysis
    Timer ptimer;
    Timer ctimer;
    Timer qtimer;

    // apply the TXN scheme on the document vectors (normalize them)
    txnScheme();

    // initialize the data arrays; keep track of the arrays locally
    ClusterData *data = new ClusterData(k, dc, wc, doc_matrix);
    float **concepts = data->concepts;
    bool *changed = data->changed;
    float *cosines = data->cosine_similarities;

    // compute initial partitioning, concepts, and quality
    initClusters(data);
    float quality = computeQ(data);
    cout << "Initial quality: " << quality << endl;


    // do spherical k-means loop
    float dQ = Q_THRESHOLD * 10;
    int iterations = 0;
    while(dQ > Q_THRESHOLD) {
        iterations++;

        // compute new clusters based on old concept vectors
        ptimer.start();
        #pragma omp parallel for
        for(int i=0; i<dc; i++) {
            int cIndx = 0;
            // only update cosine similarities if cluster has changed
            if(changed[0])
                cosines[i*k] = cosineSimilarity(data, i, 0);
            for(int j=1; j<k; j++) {
                if(changed[j]) // again, only if changed
                    cosines[i*k + j] = cosineSimilarity(data, i, j);
                if(cosines[i*k + j] > cosines[i*k + cIndx])
                    cIndx = j;
            }
            data->assignCluster(i, cIndx);
        }
        ptimer.stop();

        // update which clusters changed since last time, then swap pointers
        if(optimize)
            data->findChangedClusters();
        data->applyAssignments();

        // compute new concept vectors
        ctimer.start();
        #pragma omp parallel for
        for(int i=0; i<k; i++) {
            // only update concept vectors if cluster has changed
            if(changed[i]) {
                delete[] concepts[i];
                concepts[i] = computeConcept(data, i);
            }
        }
        ctimer.stop();

        // compute quality of new partitioning
        qtimer.start();
        float n_quality = computeQ(data);
        dQ = n_quality - quality;
        quality = n_quality;
        qtimer.stop();

        // report the quality of the current partitioning
        reportQuality(data, quality, dQ);
    }


    // report runtime statistics
    timer.stop();
    reportTime(iterations, timer.get(), ptimer.get(), ctimer.get(),
               qtimer.get());

    // return the resulting clusters and concepts in the ClusterData struct
    return data;
}
