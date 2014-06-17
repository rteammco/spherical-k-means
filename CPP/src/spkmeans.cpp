/* File: spkmeans.cpp
 *
 * Defines the basic SPKMeans functions for the generic SPKMeans class.
 */

#include "spkmeans.h"

#include <cmath>
#include <iostream>
#include <queue>
#include <vector>

#include "Galois/Timer.h"

#include "cluster_data.h"
#include "vectors.h"

using namespace std;



// Constructor: initialize variables and document norms.
SPKMeans::SPKMeans(float **doc_matrix_, int k_, int dc_, int wc_)
    : doc_matrix(doc_matrix_), k(k_), dc(dc_), wc(wc_)
{
    // init the doc vector norms
    doc_norms = new float[dc];
    for(int i=0; i<dc; i++)
        doc_norms[i] = vec_norm(doc_matrix[i], wc);

    // initialize optimization to true (optimization will happen)
    optimize = true;
}



// Destructor: clean up document norms array.
SPKMeans::~SPKMeans()
{
    delete doc_norms;
}



// Disables optimization for testing purposes.
void SPKMeans::disableOptimization()
{
    optimize = false;
}



// Enables optimization if it was disabled.
void SPKMeans::enableOptimization()
{
    optimize = true;
}



// Reports time data after running the algorithm.
void SPKMeans::reportTime(int iterations, float total_time,
                          float p_time, float c_time, float q_time)
{
    cout << "Done in " << total_time / 1000
         << " seconds after " << iterations << " iterations." << endl;
    float total = p_time + c_time + q_time;
    if(total == 0)
        cout << "No individual time stats available." << endl;
    else {
        cout << "Timers (ms): " << endl
             << "   partition [" << p_time << "] ("
                << (p_time/total)*100 << "%)" << endl
             << "   concepts [" << c_time << "] ("
                << (c_time/total)*100 << "%)" << endl
             << "   quality [" << q_time << "] ("
                << (q_time/total)*100 << "%)" << endl;
    }
}



// Applies the TXN scheme to each document vector of the given matrix.
// TXN effectively just normalizes each of the document vectors.
void SPKMeans::txnScheme()
{
    for(int i=0; i<dc; i++)
        vec_normalize(doc_matrix[i], wc);
}



// Initializes the first partition (randomly or otherwise assigned) to provide
// a starting point for the clustering algorithm.
float SPKMeans::initPartitions(ClusterData *data)
{
    // create the first arbitrary partitioning
    int split = floor(dc / k);
    cout << "Split = " << split << endl;
    int base = 1;
    for(int i=0; i<k; i++) {
        int top = base + split - 1;
        if(i == k-1)
            top = dc;

        int p_size = top - base + 1;
        data->p_sizes[i] = p_size;
        cout << "Created new partition of size " << p_size << endl;

        data->partitions[i] = new float*[p_size];
        for(int j=0; j<p_size; j++)
            data->partitions[i][j] = doc_matrix[base + j - 1];

        base = base + split;
    }

    // compute the initial concept vectors
    for(int i=0; i<k; i++) {
        data->concepts[i] =
            computeConcept(data->partitions[i], data->p_sizes[i]);
    }

    // compute initial quality of the partitions
    return computeQ(data->partitions, data->p_sizes, data->concepts);
}



// Returns the quality of the given partition by doing a dot product against
// its given concept vector.
float SPKMeans::computeQ(float **partition, int p_size, float *concept)
{
    float *sum_p = vec_sum(partition, wc, p_size);
    float quality = vec_dot(sum_p, concept, wc);
    delete sum_p;
    return quality;
}



// Returns the total quality of all partitions by summing the qualities of
// each individual partition.
float SPKMeans::computeQ(float ***partitions, int *p_sizes, float **concepts)
{
    float quality = 0;
    for(int i=0; i<k; i++)
        quality += computeQ(partitions[i], p_sizes[i], concepts[i]);
    return quality;
}



// Computes the cosine similarity value of the two given vectors (dv and cv).
// TODO - we can cache the norms of all of the document vectors
float SPKMeans::cosineSimilarity(float *cv, int doc_index)
{
    if(!optimize) {
        return vec_dot(doc_matrix[doc_index], cv, wc) /
            (vec_norm(doc_matrix[doc_index], wc) * vec_norm(cv, wc));
    }
    else {
        return vec_dot(doc_matrix[doc_index], cv, wc) /
            (doc_norms[doc_index] * vec_norm(cv, wc));
    }
}



// Computes the concept vector of the given partition. A partition is an array
// of document vectors, and the concept vector will be allocated and populated.
float* SPKMeans::computeConcept(float **partition, int p_size)
{
    float *cv = vec_sum(partition, wc, p_size);
    vec_multiply(cv, wc, (1.0 / wc));
    vec_divide(cv, wc, vec_norm(cv, wc));
    return cv;
}



// Runs the spherical k-means algorithm on the given sparse matrix D and
// clusters the data into k partitions. Non-parallel (standard) version.
ClusterData* SPKMeans::runSPKMeans()
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

    // choose an initial partitioning, and get concepts and quality
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
    while(dQ > Q_THRESHOLD) {
        iterations++;

        // compute new partitions based on old concept vectors
        ptimer.start();
        // TODO - new_partitions isn't deleted (memory leak?)
        vector<float*> *new_partitions = new vector<float*>[k];
        for(int i=0; i<k; i++)
            new_partitions[i] = vector<float*>();
        // TODO - if concept vector didn't change (i.e. partition didn't
        //        change), no need to recompute cosine similarity for
        //        any of the documents.
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
            new_partitions[cIndx].push_back(doc_matrix[i]);
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
        // TODO - if a partition didn't change since last iteration, no need
        //        to recompute its concept vector.
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
