/* File: spkmeans.cpp
 *
 * Defines the basic SPKMeans functions for the generic SPKMeans class.
 */

#include "spkmeans.h"
#include "vectors.h"



// Applies the TXN scheme to each document vector of the given matrix.
// TXN effectively just normalizes each of the document vectors.
void SPKMeans::txnScheme(float **doc_matrix, int dc, int wc)
{
    for(int i=0; i<dc; i++)
        vec_normalize(doc_matrix[i], wc);
}



// Returns the quality of the given partition by doing a dot product against
// its given concept vector.
float SPKMeans::computeQ(float **partition, int p_size, float *concept, int wc)
{
    float *sum_p = vec_sum(partition, wc, p_size);
    float quality = vec_dot(sum_p, concept, wc);
    delete sum_p;
    return quality;
}



// Returns the total quality of all partitions by summing the qualities of
// each individual partition.
float SPKMeans::computeQ(float ***partitions, int *p_sizes, float **concepts,
    int k, int wc)
{
    float quality = 0;
    for(int i=0; i<k; i++)
        quality += computeQ(partitions[i], p_sizes[i], concepts[i], wc);
    return quality;
}



// Computes the cosine similarity value of the two given vectors (dv and cv).
// TODO - we can cache the norms of all of the document vectors
float SPKMeans::cosineSimilarity(float *dv, float *cv, int wc)
{
    return vec_dot(dv, cv, wc) / (vec_norm(dv, wc) * vec_norm(cv, wc));
}



// Computes the concept vector of the given partition. A partition is an array
// of document vectors, and the concept vector will be allocated and populated.
float* SPKMeans::computeConcept(float **partition, int p_size, int wc)
{
    float *cv = vec_sum(partition, wc, p_size);
    vec_multiply(cv, wc, (1.0 / wc));
    vec_divide(cv, wc, vec_norm(cv, wc));
    return cv;
}
