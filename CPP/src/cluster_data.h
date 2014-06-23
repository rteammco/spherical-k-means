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

    // pointers to partitions and concepts
    float ***partitions;
    int *p_sizes;
    float **concepts;

    // pointers to cosine similarities, qualities, and partition change flags
    bool *changed;
    float *cValues;
    float *qualities;


    // Constructor: pass in the three required values (k, wc, dc), and set
    // partition, concept vector, and change caching pointers optionally.
    // If pointers to the lists are not provided, new lists will be
    // initialized instead.
    ClusterData(int k_, int dc_, int wc_,
            float ***ps_ = 0, int *psz_ = 0, float **cvs_ = 0,
            bool *changed_ = 0, float *cValues_ = 0, float *qualities_ = 0)
    {
        // set the variables (k, document count, word count)
        k = k_;
        dc = dc_;
        wc = wc_;

        // set partitions pointer
        if(ps_ == 0)
            partitions = new float**[k];
        else
            partitions = ps_;

        // set partition sizes pointer
        if(psz_ == 0)
            p_sizes = new int[k];
        else
            p_sizes = psz_;

        // set concepts pointer
        if(cvs_ == 0)
            concepts = new float*[k];
        else
            concepts = cvs_;

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
