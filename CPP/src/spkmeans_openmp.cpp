/* File: spkmeans_openmp.cpp
 *
 * Contains the definitions of the OpenMP version of the spherical K-means
 * algorithm. OpenMP provides multithreading.
 */

#include "spkmeans.h"

#include <cmath>
#include <iostream>
#include <mutex>
#include <queue>
#include <vector>

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
float SPKMeansOpenMP::computeQ(
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
}


// Runs the spherical k-means algorithm on the given sparse matrix D and
// clusters the data into k partitions.
ClusterData* SPKMeansOpenMP::runSPKMeans()
{
    // keep track of the run time for this algorithm
    Galois::Timer timer;
    timer.start();

    // apply the TXN scheme on the document vectors (normalize them)
    txnScheme();

    // initialize the data arrays; keep track of the arrays locally
    ClusterData *data = new ClusterData(k, dc, wc);
    float ***partitions = data->partitions;
    int *p_sizes = data->p_sizes;
    float **concepts = data->concepts;

    // compute initial quality of the partitions
    float quality = initPartitions(data);
    cout << "Initial quality: " << quality << endl;

    // keep track of all individual component times for analysis
    Galois::Timer ptimer;
    Galois::Timer ctimer;
    Galois::Timer qtimer;
    float p_time = 0;
    float c_time = 0;
    float q_time = 0;

    // do spherical k-means loop
    float dQ = Q_THRESHOLD * 10;
    int iterations = 0;
    mutex mut;
    while(dQ > Q_THRESHOLD) {
        iterations++;

        // compute new partitions based on old concept vectors
        ptimer.start();
        vector<float*> *new_partitions = new vector<float*>[k];
        #pragma omp parallel for
        for(int i=0; i<k; i++)
            new_partitions[i] = vector<float*>();
        #pragma omp parallel for
        for(int i=0; i<dc; i++) {
            int cIndx = 0;
            float cVal = cosineSimilarity(concepts[0], i);
            for(int j=1; j<k; j++) {
                float new_cVal = cosineSimilarity(concepts[j], i);
                if(new_cVal > cVal) {
                    cVal = new_cVal;
                    cIndx = j;
                }
            }
            mut.lock();
            new_partitions[cIndx].push_back(doc_matrix[i]);
            mut.unlock();
        }
        ptimer.stop();
        p_time += ptimer.get();

        // transfer the new partitions to the partitions array
        data->clearPartitions();
        for(int i=0; i<k; i++) {
            partitions[i] = new_partitions[i].data();
            p_sizes[i] = new_partitions[i].size();
        }

        // compute new concept vectors
        ctimer.start();
        data->clearConcepts();
        #pragma omp parallel for
        for(int i=0; i<k; i++)
            concepts[i] = computeConcept(partitions[i], p_sizes[i]);
        ctimer.stop();
        c_time += ctimer.get();

        // compute quality of new partitioning
        qtimer.start();
        float n_quality = computeQ(partitions, p_sizes, concepts);
        dQ = n_quality - quality;
        quality = n_quality;
        qtimer.stop();
        q_time += qtimer.get();

        cout << "Quality: " << quality << " (+" << dQ << ")" << endl;
    }


    // report runtime statistics
    timer.stop();
    reportTime(iterations, timer.get(), p_time, c_time, q_time);

    // return the resulting partitions and concepts in the ClusterData struct
    return data;
}
