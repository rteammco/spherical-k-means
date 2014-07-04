/* File: cluster_data.cpp
 *
 * Defines the ClusterData functions used to assist the SPKMeans class.
 * This class keeps track of the algorithm's data structures, and provdes
 * basic data manipulation and memory management methods.
 */

#include "cluster_data.h"



// Constructor: pass in the three required values (k, wc, dc), and setup
// the appropriate data structures.
// If pointers to the lists are not provided, new lists will be
// initialized instead.
ClusterData::ClusterData(int k_, int dc_, int wc_,
    float **concepts_, int *p_asgns_, float *doc_priorities_,
    bool *changed_, float *cValues_, float *qualities_)
{
    // set the variables (k, document count, word count)
    k = k_;
    dc = dc_;
    wc = wc_;

    // set partitions pointer, and initialize all to 0
    partitions = new float**[k];
    for(int i=0; i<k; i++)
        partitions[i] = 0;

    // set partition sizes pointer
    p_sizes = new int[k];

    // set concepts pointer
    if(concepts_== 0)
        concepts = new float*[k];
    else
        concepts = concepts_;

    // set partition assignment pointers
    if(p_asgns_ == 0)
        p_asgns = new int[dc];
    else
        p_asgns = p_asgns_;
    p_asgns_new = new int[dc];

    // set document priorities pointer
    if(doc_priorities_ == 0)
        doc_priorities = new float[dc];
    else
        doc_priorities = doc_priorities_;

    // set changed flag pointer; if newly created, initialize all to true
    if(changed_ == 0) {
        changed = new bool[k];
        for(int i=0; i<k; i++)
            changed[i] = true;
    }
    else
        changed = changed_;

    // set cosine similarity cache pointer
    if(cValues_ == 0)
        cValues = new float[k*dc];
    else
        cValues = cValues_;

    // set partition qualities cache pointer
    if(qualities_ == 0)
        qualities = new float[k];
    else
        qualities = qualities_;
}



// Destructor: calls its own clean up function
ClusterData::~ClusterData()
{
    clearMemory();
}



// Gives document doc a new partition assignment part (these are indices).
// These partitions are stored in the p_asgns_new array, and NOT the
// default one. Use swapAssignments to update.
void ClusterData::assignPartition(int doc, int partition)
{
    p_asgns_new[doc] = partition;
}



// Assigns a partition to this document (same as assignPartition above),
// but also updates the document's priority.
void ClusterData::assignPartition(int doc, int partition, float priority)
{
    assignPartition(doc, partition);
    doc_priorities[doc] = priority;
}



// Swaps the p_assignments and new_p_assignments pointers such so that
// the new values are updated without needing to manipulate memory.
void ClusterData::swapAssignments()
{
    int *temp = p_asgns;
    p_asgns = p_asgns_new;
    p_asgns_new = temp;
}



// Computes which partitions have changed, and assigns the boolean array
// accordingly.
void ClusterData::findChangedPartitions()
{
    for(int i=0; i<k; i++)
        changed[i] = false;
    for(int i=0; i<dc; i++) {
        if(p_asgns[i] != p_asgns_new[i]) {
            changed[p_asgns[i]] = true;
            changed[p_asgns_new[i]] = true;
        }
    }
}



// Cleans up partition array by deleting the list of document vectors in
// each partition slot. The vectors themselves are NOT deleted, nor is
// the actual partition list.
void ClusterData::clearPartitions()
{
    for(int i=0; i<k; i++)
        if(partitions[i] != 0)
            delete[] partitions[i];
}



// Cleans concept vector memory, permanently deleting all of the concept
// vectors. The actual concept vector list is NOT deleted.
void ClusterData::clearConcepts()
{
    for(int i=0; i<k; i++)
        delete[] concepts[i];
}



// Cleans up all structs that can be cleaned, and frees up all used memory.
// WARNING: this will turn all data structure pointers to NULL.
void ClusterData::clearMemory()
{
    // clean up partition, size, and concept vector arrays
    if(partitions != 0) {
        clearPartitions();
        delete[] partitions;
        partitions = 0;
    }
    if(p_sizes != 0) {
        delete[] p_sizes;
        p_sizes = 0;
    }

    // clean up concept vectors
    if(concepts != 0) {
        clearConcepts();
        delete[] concepts;
        concepts = 0;
    }

    // clean up partition assignment arrays
    if(p_asgns != 0)
        delete[] p_asgns;
    if(p_asgns_new != 0)
        delete[] p_asgns_new;

    // clean up document priorities array
    if(doc_priorities != 0)
        delete[] doc_priorities;

    // clean up change cache arrays
    if(changed != 0) {
        delete[] changed;
        changed = 0;
    }
    if(cValues != 0) {
        delete[] cValues;
        cValues = 0;
    }
    if(qualities != 0) {
        delete[] qualities;
        qualities = 0;
    }
}

