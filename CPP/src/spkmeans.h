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
    // matrix setup schemes
    void txnScheme(float **doc_matrix, int dc, int wc);

    // quality computation functions
    float computeQ(float **partition, int p_size, float *concept, int wc);
    virtual float computeQ(
        float ***partitions, int *p_sizes, float **concepts, int k, int wc);

    // spkmeans algorithm computation function
    float cosineSimilarity(float *dv, float *cv, int wc);
    float* computeConcept(float **partition, int p_size, int wc);

  public:
    // the algorithm is implemented differently by each type of paradigm
    virtual ClusterData* runSPKMeans(float **doc_matrix, int k, int dc, int wc);
};



// OpenMP version of the SPKMeans algorithm
class SPKMeansOpenMP : SPKMeans {
  private:
    unsigned int num_threads;
    float computeQ(float ***partitions, int *p_sizes, float **concepts,
        int k, int wc);

  public:
    // constructor: set the number of threads
    SPKMeansOpenMP(unsigned int t_ = 1);

    // returns the actual number of threads Galois will use
    unsigned int getNumThreads();

    // run the algorithm
    ClusterData* runSPKMeans(float **doc_matrix, int k, int dc, int wc);
};



// Galois version of the SPKMeans algorithm
class SPKMeansGalois : SPKMeans {
  private:
    unsigned int num_threads;

  public:
    // constructor: set the number of threads and initialize Galois
    SPKMeansGalois(unsigned int t_ = 1);

    // returns the actual number of threads Galois will use
    unsigned int getNumThreads();

    // run the algorithm
    ClusterData* runSPKMeans(float **doc_matrix, int k, int dc, int wc);
};



#endif
