/* File: cluster_data.h
 *
 * Contains a struct for the clustering algorithm(s) to use to easily store
 * and return the partitions, partition sizes, and concept vectors, along
 * with other variables.
 */

#ifndef CLUSTER_DATA_H
#define CLUSTER_DATA_H


// ClusterData struct can contain partition and concept vector pointers, and
// functions to clean out the memory.
struct ClusterData
{
    // clustering variables (k, word count, document count)
    int k;
    int dc;
    int wc;

    // pointers to partitions (OLD CODE)
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


    // Constructor: pass in the three required values (k, wc, dc), and set
    // partition, concept vector, and change caching pointers optionally.
    // If pointers to the lists are not provided, new lists will be
    // initialized instead.
    ClusterData(int k_, int dc_, int wc_,
            float **cvs_ = 0, int *p_asgns_ = 0, float *doc_priorities_ = 0,
            bool *changed_ = 0, float *cValues_ = 0, float *qualities_ = 0)
    {
        // set the variables (k, document count, word count)
        k = k_;
        dc = dc_;
        wc = wc_;

        // set partitions pointer
        partitions = new float**[k];

        // set partition sizes pointer
        p_sizes = new int[k];

        // set concepts pointer
        if(cvs_ == 0)
            concepts = new float*[k];
        else
            concepts = cvs_;

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
    ~ClusterData()
    {
        clearMemory();
    }


    // Gives document doc a new partition assignment part (these are indices).
    // These partitions are stored in the p_asgns_new array, and NOT the
    // default one. Use swapAssignments to update.
    // This function also automatically updates the partition change values.
    void assignPartition(int doc, int partition)
    {
        p_asgns_new[doc] = partition;
        if(p_asgns[doc] != partition) {
            changed[partition] = true;
            changed[p_asgns[doc]] = true;
        }
    }


    // Assigns a partition to this document (same as assignPartition above),
    // but also updates the document's priority.
    void assignPartition(int doc, int partition, float priority)
    {
        assignPartition(doc, partition);
        doc_priorities[doc] = priority;
    }


    // Swaps the p_assignments and new_p_assignments pointers such so that
    // the new values are updated without needing to manipulate memory.
    void swapAssignments()
    {
        int *temp = p_asgns;
        p_asgns = p_asgns_new;
        p_asgns_new = temp;
    }


    // Cleans up partition array by deleting the list of document vectors in
    // each partition slot. The vectors themselves are NOT deleted, nor is
    // the actual partition list.
    void clearPartitions()
    {
        for(int i=0; i<k; i++)
            delete[] partitions[i];
    }


    // Cleans concept vector memory, permanently deleting all of the concept
    // vectors. The actual concept vector list is NOT deleted.
    void clearConcepts()
    {
        for(int i=0; i<k; i++)
            delete[] concepts[i];
    }


    // clean up the partitions and concept vector pointers
    void clearMemory()
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
        if(concepts != 0) {
            clearConcepts();
            delete[] concepts;
            concepts = 0;
        }

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

};


#endif
