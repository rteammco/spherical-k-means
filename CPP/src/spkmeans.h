/* File: spkmeans.h
 *
 * Contains the SPKMeans class (and generalizations to different parallel
 * paradigms). The SPKMeans classes do the actual spherical k-means algorithm
 * computation and contain helper functions for the algorithm.
 */

#ifndef SPKMEANS_H
#define SPKMEANS_H

#include "cluster_data.h"

#define Q_THRESHOLD 0.001



// Abstract implementation of the SPKMeans algorithm
class SPKMeans {
  protected:
    // clustering variables
    float **doc_matrix;
    int k;
    int dc;
    int wc;

    // cache for document norms so they're not recomputed each time
    float *doc_norms;

    // optimization flag
    bool optimize;
    
    // report timer stats
    void reportTime(int iterations, float total_time,
                    float p_time = 0, float c_time = 0, float q_time = 0);

    // matrix setup schemes
    void txnScheme();

    // initial partition setup
    void initPartitions(ClusterData *data);

    // quality computation functions
    float computeQ(float **partition, int p_size, float *concept);
    float computeQ(ClusterData *data);

    // spkmeans algorithm computation function
    float* computeConcept(float **partition, int p_size);

  public:
    // initialize wc, dc, k, and doc_matrix, and document norms
    SPKMeans(float **doc_matrix_, int k_, int dc_, int wc_);
    // clean up memory
    ~SPKMeans();

    // switches for optimization
    void disableOptimization();
    void enableOptimization();

    // cosine similarity compute function: made public for use by Galois
    // computation structs
    float cosineSimilarity(float *cv, int doc_index);

    // the algorithm is implemented differently by each type of paradigm
    virtual ClusterData* runSPKMeans();
};



// OpenMP version of the SPKMeans algorithm
class SPKMeansOpenMP : public SPKMeans {
  private:
    unsigned int num_threads;
    //float computeQ(float ***partitions, int *p_sizes, float **concepts);

  public:
    // constructor: set the number of threads
    SPKMeansOpenMP(float **doc_matrix_, int k_, int dc_, int wc_,
        unsigned int t_ = 1);

    // returns the actual number of threads Galois will use
    unsigned int getNumThreads();

    // run the algorithm
    ClusterData* runSPKMeans();
};



// Galois version of the SPKMeans algorithm
class SPKMeansGalois : public SPKMeans {
  private:
    unsigned int num_threads;

  public:
    // constructor: set the number of threads and initialize Galois
    SPKMeansGalois(float **doc_matrix_, int k_, int dc_, int wc_,
        unsigned int t_ = 1);

    // returns the actual number of threads Galois will use
    unsigned int getNumThreads();

    // run the algorithm
    ClusterData* runSPKMeans();
};



#endif
