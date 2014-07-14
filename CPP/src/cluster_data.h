/* File: cluster_data.h
 *
 * Contains a class for the clustering algorithm(s) to use and easily store
 * and return the partition assignments, caching arrays, concept vectors,
 * and other data structures used by spherical k-means.
 */

#ifndef CLUSTER_DATA_H
#define CLUSTER_DATA_H

#include <vector>


// This struct is used to store a word value and index pair, used by the
// Document struct to map words.
struct ValueIndexPair {
    float value;
    int index;
};


// This struct maps each document to the value-index pairs of all words
// associated with it.
struct Document {
    int count;
    std::vector<ValueIndexPair> words;
};


// ClusterData class can contain partition assignments and concept vector
// pointers, and functions to manage optimizations and memory.
class ClusterData {

  public:

    // clustering variables (k, word count, document count)
    int k;
    int dc;
    int wc;

    // pointers to concept vectors and partition assignments (current and new)
    int *p_asgns;
    int *p_asgns_new;
    float **concepts;

    // document priority and heuristic tracking variables
    float *doc_priorities;
    float total_priority;
    float total_moved_priority;
    int num_moved;

    // pointers to cosine similarities, qualities, and cluster change flags
    bool *changed;
    float *cosine_similarities;
    float *qualities;

    // document data structures that map documents to words
    Document *docs;


    // Constructor: sets up variables and data structures.
    ClusterData(int k_, int dc_, int wc_, float **doc_matrix,
            float **cvs_ = 0, int *p_asgns_ = 0, float *doc_priorities_ = 0,
            bool *changed_ = 0, float *cosine_similarities_ = 0,
            float *qualities_ = 0);

    // Destructor: calls its own clean up function.
    ~ClusterData();

    // Assigns a cluster to the given document.
    void assignCluster(int doc, int cluster);

    // Assigns a cluster and priority to the given document.
    void assignCluster(int doc, int cluster, float priority);

    // Swaps new assignments for the default ones (updates the assignments).
    void applyAssignments();

    // Returns the average priority of all documents.
    float getAveragePriority();

    // Returns the average priority of all documents that have moved.
    float getAverageMovedPriority();

    // Returns the average priority of all documents that have NOT moved.
    float getAverageStayPriority();

    // Updates which clusters have been changed since last partitioning.
    void findChangedClusters();

    // Cleans concept vector memory.
    void clearConcepts();

    // Clean up all data memory.
    void clearMemory();

};


#endif
