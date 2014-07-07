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
        data->concepts[i] = temp_computeConcept(data, i);
}



// Returns the quality of the given partition by doing a dot product against
// its given concept vector.
float SPKMeans::computeQ(float **partition, int p_size, float *concept)
{
    float *sum_p = vec_sum(partition, wc, p_size);
    float quality = vec_dot(sum_p, concept, wc);
    delete[] sum_p;
    return quality;
}



// Returns the total quality of all partitions by summing the qualities of
// each individual partition. If optimization is enabled, uses cached values
// whenever possible.
float SPKMeans::computeQ(ClusterData *data)
{
    float quality = 0;

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
float SPKMeans::cosineSimilarity(float *cv, int doc_index)
{
    if(!optimize) { // if optimization is off, don't use caching
        return vec_dot(doc_matrix[doc_index], cv, wc) /
            (vec_norm(doc_matrix[doc_index], wc) * vec_norm(cv, wc));
    }
    else {
        return vec_dot(doc_matrix[doc_index], cv, wc) /
            (doc_norms[doc_index] * vec_norm(cv, wc));
    }
}



// Determines which partitions have changed between the old partitioning
// (in the ClusterData array) and the new partitions (in the given vector).
// Sets the ClusterData::changed values acoordingly.
// For parallel implementations, this needs to be adjusted to be order-
// invariant. This code assumes documents were inserted into partitions in
// order of their indices.
void SPKMeans::findChangedPartitions(vector<float*> *new_partitions,
    ClusterData *data)
{
    if(optimize) { // skip if not optimizing
        for(int i=0; i<(data->k); i++) {
            if(data->p_sizes[i] == new_partitions[i].size()) {
                data->changed[i] = false;
                // only check forward direction: documents were inserted
                // in order by index
                for(int j=0; j<(data->p_sizes[i]); j++) {
                    if(data->partitions[i][j] != new_partitions[i][j]) {
                        data->changed[i] = true;
                        break;
                    }
                }
            }
            else
                data->changed[i] = true;
        }
    }
}



// Same as findChangedPartitions, but without the assumption that the new
// partitions were computed in order. Use this version when computing the
// partitions in parallel or when unsure what the order of computation was.
void SPKMeans::findChangedPartitionsUnordered(vector<float*> *new_partitions,
    ClusterData *data)
{
    if(optimize) { // skip if not optimizing
        for(int i=0; i<(data->k); i++) {
            // for each partition, check if new and old are same size
            if(data->p_sizes[i] == new_partitions[i].size()) {
                data->changed[i] = false;
                // for each document in old partition
                for(int j=0; j<(data->p_sizes[i]); j++) {
                    // check if document is in new partition (we need
                    // std::find here because order is not guaranteed due
                    // to the parallelization)
                    if(find(new_partitions[i].begin(),
                            new_partitions[i].end(),
                            data->partitions[i][j]) == new_partitions[i].end())
                    {
                        data->changed[i] = true;
                        break;
                    }
                }
            }
            else
                data->changed[i] = true;
        }
    }
}



// Clears the partitions in the ClusterData struct, and copies the new
// partitions (and size information) over from the given vector.
void SPKMeans::copyPartitions(vector<float*> *new_partitions,
    ClusterData *data)
{
    data->clearPartitions();
    for(int i=0; i<(data->k); i++) {
        data->p_sizes[i] = new_partitions[i].size();
        data->partitions[i] = new_partitions[i].data();
    }
}



// Computes the concept vector of the given partition. A partition is an array
// of document vectors, and the concept vector will be allocated and populated.
float* SPKMeans::computeConcept(float **partition, int p_size)
{
    float *cv = vec_sum(partition, wc, p_size);
    vec_divide(cv, wc, wc);
    vec_normalize(cv, wc);
    return cv;
}



/***** TEMP ADDED */
float* SPKMeans::temp_computeConcept(ClusterData *data, int pIndx)
{
    float *concept = new float[wc];
    for(int i=0; i<wc; i++)
        concept[i] = 0;
    for(int i=0; i<dc; i++) { // for each document
        if(data->p_asgns[i] == pIndx) { // if doc in cluster
            for(int j=0; j<wc; j++) // add words to concept
                concept[j] += doc_matrix[i][j];
        }
    }
    vec_divide(concept, wc, wc);
    vec_normalize(concept, wc);
    return concept;
}

// Runs the spherical k-means algorithm on the given sparse matrix D and
// clusters the data into k partitions. Non-parallel (standard) version.
ClusterData* SPKMeans::runSPKMeans()
{
    // TODO - temporary (implement this in reader)
    // convert document data format!
    txnScheme();
    temp_Document* docs[dc];
    for(int i=0; i<dc; i++) {
        vector<temp_ValueIndex*> non_zeros;
        for(int j=0; j<wc; j++) {
            if(doc_matrix[i][j] > 0) {
                temp_ValueIndex *vi = new temp_ValueIndex();
                vi->value = doc_matrix[i][j];
                vi->index = j;
                non_zeros.push_back(vi);
            }
        }
        docs[i] = new temp_Document();
        int num_nz = non_zeros.size();
        docs[i]->num_nonzero = num_nz;
        docs[i]->non_zeros = new temp_ValueIndex*[num_nz];
        for(int j=0; j<num_nz; j++) {
            temp_ValueIndex *vi = non_zeros[j];
            docs[i]->non_zeros[j] = non_zeros[j];
        }
    }

    // keep track of the run time for this algorithm
    Galois::Timer timer;
    timer.start();

    // apply the TXN scheme on the document vectors (normalize them)
    txnScheme();

    // initialize the data arrays; keep track of the arrays locally
    ClusterData *data = new ClusterData(k, dc, wc);
    float **concepts = data->concepts;
    bool *changed = data->changed;
    float *cValues = data->cValues;
    float *qualities = data->qualities;

    // compute initial partitions, concepts, and quality
    initPartitions(data);
    float quality = computeQ(data);
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
        for(int i=0; i<dc; i++) {
            /*// TODO - temp added
            int pIndx = 0;
            float temp = 0;
            for(int j=0; j<k; j++) {
                if(changed[0]) {
                    float cnorm = vec_norm(concepts[j], wc); // TODO - same op. for cvs
                    float dnorm = doc_norms[i]; // TODO - we can cache this better too
                    // here is where we save time: compute the dot product!
                    float dotp = 0;
                    for(int a=0; a<(docs[i]->num_nonzero); a++) {
                        int word = docs[i]->non_zeros[a]->index;
                        float value = docs[i]->non_zeros[a]->value;
                        dotp += concepts[j][word] * value;
                    }
                    cValues[i*k + j] = dotp / (cnorm * dnorm);
                }
                if(cValues[i*k + j] > cValues[i*k + pIndx]) {
                    pIndx = j;
                    temp = cValues[i*k + j];
                }
            }*/

            // only update cosine similarities if partitions have changed
            // or if optimization is disabled
            int pIndx2 = 0;
            if(changed[0])
                cValues[i*k] = cosineSimilarity(concepts[0], i);
            for(int j=1; j<k; j++) {
                if(changed[j]) // again, only if changed
                    cValues[i*k + j] = cosineSimilarity(concepts[j], i);
                if(cValues[i*k + j] > cValues[i*k + pIndx2])
                    pIndx2 = j;
            }
            //if(pIndx != pIndx2)
            //    cout << "ERROR: FAIL: " << cValues[i*k + pIndx2] << " vs. " << temp << endl;
            data->assignPartition(i, pIndx2);
        }
        ptimer.stop();
        p_time += ptimer.get();

        // update which partitions changed since last time, then swap pointers
        data->findChangedPartitions();
        data->swapAssignments();

        // compute new concept vectors
        ctimer.start();
        for(int i=0; i<k; i++) {
            if(changed[i]) {
                delete[] concepts[i];
                concepts[i] = temp_computeConcept(data, i);
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

    delete data;
    // return the resulting partitions and concepts in the ClusterData struct
    return 0;
}
