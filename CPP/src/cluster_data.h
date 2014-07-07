/* File: cluster_data.h
 *
 * Contains a class for the clustering algorithm(s) to use and easily store
 * and return the partition assignments, caching arrays, concept vectors,
 * and other data structures used by spherical k-means.
 */

#ifndef CLUSTER_DATA_H
#define CLUSTER_DATA_H


// ClusterData class can contain partition and concept vector pointers, and
// functions to clean out the memory.
class ClusterData {

  public:

    // clustering variables (k, word count, document count)
    int k;
    int dc;
    int wc;

    // TODO - remove this
    float ***partitions;
    int *p_sizes;

    // pointers to concept vectors and partition assignments (current and new),
    // and document priorities
    int *p_asgns;
    int *p_asgns_new;
    float **concepts;
    float *doc_priorities;

    // pointers to cosine similarities, qualities, and partition change flags
    bool *changed;
    float *cValues;
    float *qualities;


    // Constructor: sets up variables and data structures.
    ClusterData(int k_, int dc_, int wc_,
            float **cvs_ = 0, int *p_asgns_ = 0, float *doc_priorities_ = 0,
            bool *changed_ = 0, float *cValues_ = 0, float *qualities_ = 0);

    // Destructor: calls its own clean up function.
    ~ClusterData();

    // Assigns a partition to the given document.
    void assignPartition(int doc, int partition);

    // Assigns a partition and priority to the given document.
    void assignPartition(int doc, int partition, float priority);

    // Swaps new assignments for the default ones (updates the assignments).
    void swapAssignments();

    // Updates which partitions have been changed since last partitioning.
    void findChangedPartitions();

    // Cleans up partition array. TODO - remove
    void clearPartitions();

    // Cleans concept vector memory.
    void clearConcepts();

    // Clean up all data memory.
    void clearMemory();

};


#endif
