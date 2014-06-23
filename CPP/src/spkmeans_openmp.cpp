/* File: spkmeans_openmp.cpp
 *
 * Contains the definitions of the OpenMP version of the spherical K-means
 * algorithm. OpenMP provides multithreading.
 */

#include "spkmeans.h"

#include <algorithm>
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

    // apply the TXN scheme on the document vectors (normalize them)
    txnScheme();

    // initialize the data arrays; keep track of the arrays locally
    ClusterData *data = new ClusterData(k, dc, wc);
    float ***partitions = data->partitions;
    int *p_sizes = data->p_sizes;
    float **concepts = data->concepts;

    // choose an initial partitioning, and get first concepts
    initPartitions(data);

    // keep track of all individual component times for analysis
    Galois::Timer ptimer;
    Galois::Timer ctimer;
    Galois::Timer qtimer;
    float p_time = 0;
    float c_time = 0;
    float q_time = 0;

    // keep track of which partitions were changed, the cosine similarities
    // between all docs and clusters, and the qualities of each cluster.
    bool changed[k];
    for(int i=0; i<k; i++)
        changed[i] = true;
    float cValues[k*dc];
    float qualities[k];

    // compute initial quality, and cache the quality values
    float quality = computeQ(data);
    cout << "Initial quality: " << quality << endl;

    // do spherical k-means loop
    float dQ = Q_THRESHOLD * 10;
    int iterations = 0;
    mutex mut;
    while(dQ > Q_THRESHOLD) {
        iterations++;

        // compute new partitions based on old concept vectors
        ptimer.start();
        // TODO - new_partitions isn't deleted (memory leak?)
        vector<float*> *new_partitions = new vector<float*>[k];
        #pragma omp parallel for
        for(int i=0; i<k; i++)
            new_partitions[i] = vector<float*>();
        #pragma omp parallel for
        for(int i=0; i<dc; i++) {
            int cIndx = 0;
            // only update cosine similarities if partitions have changed
            // or if optimization is disabled
            if(changed[0] || !optimize)
                cValues[i*k] = cosineSimilarity(concepts[0], i);
            for(int j=1; j<k; j++) {
                if(changed[j] || !optimize) // again, only if changed
                    cValues[i*k + j] = cosineSimilarity(concepts[j], i);
                if(cValues[i*k + j] > cValues[i*k + cIndx])
                    cIndx = j;
            }
            mut.lock();
            new_partitions[cIndx].push_back(doc_matrix[i]);
            mut.unlock();
        }
        ptimer.stop();
        p_time += ptimer.get();

        // check if partitions changed since last time
        for(int i=0; i<k; i++) {
            // for each partition, check if new and old are same size
            if(p_sizes[i] == new_partitions[i].size()) {
                changed[i] = false;
                // for each document in old partition
                for(int j=0; j<p_sizes[i]; j++) {
                    // check if document is in new partition (we need std::find
                    // here because order is not guaranteed due to the
                    // parallelization.
                    if(find(new_partitions[i].begin(),
                            new_partitions[i].end(),
                            partitions[i][j]) == new_partitions[i].end())
                    {
                        changed[i] = true;
                        break;
                    }
                }
            }
            else
                changed[i] = true;
        }

        // transfer the new partitions to the partitions array
        data->clearPartitions();
        for(int i=0; i<k; i++) {
            partitions[i] = new_partitions[i].data();
            p_sizes[i] = new_partitions[i].size();
        }

        // compute new concept vectors
        ctimer.start();
        #pragma omp parallel for
        for(int i=0; i<k; i++) {
            // only update concept vectors if partition has changed
            if(changed[i] || !optimize) {
                delete[] concepts[i];
                concepts[i] = computeConcept(partitions[i], p_sizes[i]);
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

        cout << "Quality: " << quality << " (+" << dQ << ")";// << endl;
        // TODO - TEMP debug message
        int num_same = 0;
        for (int i=0; i<k; i++)
            if(!changed[i])
                num_same++;
        cout << " --- " << num_same << " partitions are the same." << endl;
    }


    // report runtime statistics
    timer.stop();
    reportTime(iterations, timer.get(), p_time, c_time, q_time);

    // return the resulting partitions and concepts in the ClusterData struct
    return data;
}
