/* File: spkmeans.cpp
 *
 * Defines the basic SPKMeans functions for the generic SPKMeans class.
 * This class provides a single-threaded implementation of the algorithm,
 * as well as some basic functions useful for the multithreaded versions.
 */

#include "spkmeans.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector> // TODO - debug

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

    // default scheme to TXN
    prep_scheme = TXN_SCHEME;

    // initialize optimization to true (optimization will happen)
    optimize = true;
}



// Destructor: clean up document norms array.
SPKMeans::~SPKMeans()
{
    delete[] doc_norms;
}



// Set the current scheme to the given type.
void SPKMeans::setScheme(SPKMeans::Scheme type)
{
    prep_scheme = type;
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



// Reports the overall quality and, if optimizing, also displays how many
// partitions have changed.
void SPKMeans::reportQuality(ClusterData *data, float quality, float dQ)
{
    cout << "Quality: " << quality << " (+" << dQ << ")";
    if(optimize) {
        int num_same = 0;
        for (int i=0; i<k; i++)
            if(!(data->changed[i]))
                num_same++;
        cout << " --- " << num_same << " partitions are the same." << endl;
    }
    else
        cout << " --- optimization disabled." << endl;
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
// TXN scheme must be set, otherwise this function will do nothing.
void SPKMeans::txnScheme()
{
    if(prep_scheme != SPKMeans::TXN_SCHEME)
        return;

    for(int i=0; i<dc; i++)
        vec_normalize(doc_matrix[i], wc);
}



// Initializes the first partition (randomly or otherwise assigned) to provide
// a starting point for the clustering algorithm.
void SPKMeans::initPartitions(ClusterData *data)
{
    // choose an initial partitioning
    int split = floor(dc / k);
    cout << "Split = " << split << endl;
    int base = 1;
    for(int i=0; i<k; i++) {
        int top = base + split - 1;
        if(i == k-1)
            top = dc;
        int p_size = top - base + 1;
        for(int j=0; j<p_size; j++)
            data->p_asgns[base-1 + j] = i;
        base = base + split;
    }

    // compute the initial concept vectors
    for(int i=0; i<k; i++)
        data->concepts[i] = computeConcept(data, i);
}



// Returns the total quality of all partitions by summing the qualities of
// each individual partition. If optimization is enabled, uses cached values
// whenever possible.
float SPKMeans::computeQ(ClusterData *data)
{
    float quality = 0;

    // NOTE - there was a single-partition computeQ function,
    //        re-implement or not?
    for(int i=0; i<k; i++) {
        if(data->changed[i]) {
            // initialize sum vector to 0
            float sum_p[wc];
            for(int a=0; a<wc; a++)
                sum_p[a] = 0;
            // add all documents associated with this partition
            for(int j=0; j<dc; j++) {
                if(data->p_asgns[j] == i) {
                    for(int a=0; a<wc; a++)
                        sum_p[a] += doc_matrix[j][a];
                }
            }
            data->qualities[i] = vec_dot(sum_p, data->concepts[i], wc);
        }
        quality += data->qualities[i];
    }

    return quality;
}



// Computes the cosine similarity value of the two given vectors (dv and cv).
float SPKMeans::cosineSimilarity(ClusterData *data, int doc_index, int cluster)
{
    // TODO - same optimization for concepts? can we cache doc norms better?
    float cnorm = vec_norm(data->concepts[cluster], wc); // segfault??
    float dnorm = doc_norms[doc_index];

    // here is where we save time: compute the dot product!
    float dotp = 0;
    for(int i=0; i<(data->docs[doc_index].count); i++) {
        int word = data->docs[doc_index].words[i].index;
        float value = data->docs[doc_index].words[i].value;
        dotp += data->concepts[cluster][word] * value;
    }
    
    return dotp / (dnorm * cnorm);
}



// Computes the concept vector of the given partition (by index). The
// partition documents are accessed using the ClusterData struct, and the
// associated concept vector will be allocated and populated.
float* SPKMeans::computeConcept(ClusterData *data, int pIndx)
{
    // create the concept vector and initialize it to 0
    float *concept = new float[wc];
    for(int i=0; i<wc; i++)
        concept[i] = 0;

    // find documents associated w/ this partition, and sum them
    for(int i=0; i<dc; i++) { // for each document
        if(data->p_asgns[i] == pIndx) { // if doc in cluster
            for(int j=0; j<wc; j++) // add words to concept
                concept[j] += doc_matrix[i][j];
        }
    }

    // normalize the concept vector and return it
    vec_divide(concept, wc, wc);
    vec_normalize(concept, wc);
    return concept;
}



// Runs the spherical k-means algorithm on the given sparse matrix D and
// clusters the data into k partitions. Non-parallel (standard) version.
ClusterData* SPKMeans::runSPKMeans()
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
        for(int i=0; i<dc; i++) {
            // only update cosine similarities if partitions have changed
            // or if optimization is disabled
            int pIndx = 0;
            if(changed[0])
                cValues[i*k] = cosineSimilarity(data, i, 0);
            for(int j=1; j<k; j++) {
                if(changed[j]) // again, only if changed
                    cValues[i*k + j] = cosineSimilarity(data, i, j);
                if(cValues[i*k + j] > cValues[i*k + pIndx])
                    pIndx = j;
            }
            data->assignPartition(i, pIndx);
        }
        ptimer.stop();
        p_time += ptimer.get();

        // update which partitions changed since last time, then swap pointers
        if(optimize)
            data->findChangedPartitions();
        data->swapAssignments();

        // compute new concept vectors
        ctimer.start();
        for(int i=0; i<k; i++) {
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
