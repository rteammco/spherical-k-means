/* File: cluster_data.cpp
 *
 * Defines the ClusterData functions used to assist the SPKMeans class.
 * This class keeps track of the algorithm's data structures, and provdes
 * basic data manipulation and memory management methods.
 */

#include "cluster_data.h"

#include <vector>



// Constructor: pass in the four required values (k, wc, dc, doc_matrix), and
// set up the appropriate data structures.
// If pointers to the optional lists are not provided, new lists will be
// initialized instead.
ClusterData::ClusterData(int k_, int dc_, int wc_, float **doc_matrix,
    float **concepts_, int *p_asgns_, float *doc_priorities_,
    bool *changed_, float *cosine_similarities_, float *qualities_)
{
    // set the size variables (k, document count, word count)
    k = k_;
    dc = dc_;
    wc = wc_;

    // initialize and fill document data structure
    docs = new Document[dc];
    for(int i=0; i<dc; i++) {
        // find all non-zero words in this document
        Document d;
        for(int j=0; j<wc; j++) {
            if(doc_matrix[i][j] > 0) {
                ValueIndexPair vi;
                vi.value = doc_matrix[i][j];
                vi.index = j;
                d.words.push_back(vi);
            }
        }
        d.count = d.words.size();
        docs[i] = d;
    }

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
    if(cosine_similarities_ == 0)
        cosine_similarities = new float[k*dc];
    else
        cosine_similarities = cosine_similarities_;

    // set cluster qualities cache pointer
    if(qualities_ == 0)
        qualities = new float[k];
    else
        qualities = qualities_;

    // init all counters to 0
    total_priority = 0;
    total_moved_priority = 0;
    num_moved = 0;
}



// Destructor: calls its own clean up function
ClusterData::~ClusterData()
{
    clearMemory();
}



// Gives document doc a new cluster assignment part (these are indices).
// These clusters are stored in the p_asgns_new array, and NOT the
// default one. Use swapAssignments to update.
void ClusterData::assignCluster(int doc, int cluster)
{
    p_asgns_new[doc] = cluster;
}



// Assigns a cluster to this document (same as assignCluster above),
// but also updates the document's priority and moved statistics.
void ClusterData::assignCluster(int doc, int cluster, float priority)
{
    p_asgns_new[doc] = cluster;
    doc_priorities[doc] = priority;
    total_priority += priority;

    // if new assignment is different, this document moved
    if(cluster != p_asgns[doc]) {
        num_moved++;
        total_moved_priority += priority;
    }
}



// Swaps the p_assignments and new_p_assignments pointers such so that
// the new values are updated without needing to manipulate memory.
// This also resets the priority and moved counts to 0.
void ClusterData::applyAssignments()
{
    int *temp = p_asgns;
    p_asgns = p_asgns_new;
    p_asgns_new = temp;
    total_priority = 0;
    total_moved_priority = 0;
    num_moved = 0;
}



// Returns the average priority of all documents.
float ClusterData::getAveragePriority()
{
    return total_priority / dc;
}



// Returns the average priority of all documents that were moved.
float ClusterData::getAverageMovedPriority()
{
    if(num_moved > 0)
        return total_moved_priority / num_moved;
    else
        return 0;
}



// Returns the average priority of all documents that were NOT moved.
float ClusterData::getAverageStayPriority()
{
    if(num_moved < dc)
        return (total_priority - total_moved_priority) / (dc - num_moved);
    else
        return 0;
}



// Computes which clusters have changed, and assigns the boolean array
// accordingly.
void ClusterData::findChangedClusters()
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
    // clean up the document data structures
    if(docs != 0) {
        delete[] docs;
        docs = 0;
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
    if(cosine_similarities != 0) {
        delete[] cosine_similarities;
        cosine_similarities = 0;
    }
    if(qualities != 0) {
        delete[] qualities;
        qualities = 0;
    }
}
