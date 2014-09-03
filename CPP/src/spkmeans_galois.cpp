/* File: spkmeans_galois.cpp
 *
 * Defines the Galois version of the SKPMeans class. This version directly
 * uses the Galois parallel library for multithreaded computation to speed
 * up the algorithm's runtime on multicore systems.
 */

#include "spkmeans.h"

#include "timer.h"

#include <functional>
#include <iostream>

#include "Galois/Galois.h"
#include "Galois/Graph/Graph.h"
#include "llvm/ADT/SmallVector.h"

#include <boost/iterator/counting_iterator.hpp>

using namespace std;



// Constructor: set number of threads and initialize Galois.
SPKMeansGalois::SPKMeansGalois(
    float **doc_matrix_, int k_, int dc_, int wc_, unsigned int t_)
    : SPKMeans::SPKMeans(doc_matrix_, k_, dc_, wc_)
{
    // if number of threads given is <= 0, set to max
    if(t_ <= 0)
        t_ = sysconf(_SC_NPROCESSORS_ONLN);
        // TODO - no support in Galois API for this?
    Galois::setActiveThreads(t_);
    num_threads = Galois::getActiveThreads();
}


// Returns the actual number of threads that Galois will use.
unsigned int SPKMeansGalois::getNumThreads()
{
    return num_threads;
}



/*// Scrap code (for the potential Graph implementaton)
struct DataNode {
    int id;
};
struct WordNode : public DataNode {
    float value;
};
typedef Galois::Graph::LC_CSR_Graph<DataNode, float> Graph;
struct SPKMeansRun {

    // copy this stuff from the current ClusterData
    ClusterData *data;
    bool *changed;
    float *cosines;
    float *has_docs;
    int k;

    void operator() (Graph::GraphNode doc,
                     Galois::UserContext<Graph::GraphNode>& ctx)
    {
        int i = 0;//doc->index;
        int cIndx = 0;
        if(changed[0])
            cosines[i*k] = 0;//cosineSimilarity(data, i, 0);
        for(int j=0; j<k; j++) {
            if(changed[j])
                cosines[i*k + j] = 0;//cosineSimilarity(data, i, j);
            if(cosines[i*k + j] > cosines[i*k + cIndx])
                cIndx = j;
        }
        data->assignCluster(i, cIndx);//, priority);
        has_docs[cIndx] = true;

        // TODO - now we have to do online recompuation for the concept vector
        // and quality... this will require a thread lock :(
    }

};*/



// Runs SPKMeans in the same parallel manner as OpenMP, without using any
// additional online schemes or priorities.
struct ComputeClustersBasic {

    ClusterData *data;

    function<float(ClusterData*, int, int)> cosineSimilarity;

    // Constructor: assign the ClusterData pointer
    ComputeClustersBasic(ClusterData *data_) : data(data_) { }

    // Galois operator: run the clustering computation
    void operator() (int &i, Galois::UserContext<int> &ctx)
    {
        int k = data->k;
        bool *changed = data->changed;
        float *cosines = data->cosine_similarities;

        // find the cluster with the best cosine similarity, and assign it
        int cIndx = 0;
        if(changed[0])
            cosines[i*k] = cosineSimilarity(data, i, 0);
        for(int j=1; j<k; j++) {
            if(changed[j])
                cosines[i*k + j] = cosineSimilarity(data, i, j);
            if(cosines[i*k + j] > cosines[i*k + cIndx])
                cIndx = j;
        }
        data->assignCluster(i, cIndx);
    }
};



// Runs SPKMeans in with the online algorithm to loop constantly until
// convergence. This version can make use of the priority function.
struct ComputeClustersOnline : public ComputeClustersBasic {

    // Galois operator: run the clustering computation, with online additions
    void operator() (int &i, Galois::UserContext<int> &ctx)
    {
        int k = data->k;
        bool *changed = data->changed;
        float *cosines = data->cosine_similarities;

        // find the cluster with the best cosine similarity, and assign it
        int cIndx = 0;
        if(changed[0])
            cosines[i*k] = cosineSimilarity(data, i, 0);
        for(int j=1; j<k; j++) {
            if(changed[j])
                cosines[i*k + j] = cosineSimilarity(data, i, j);
            if(cosines[i*k + j] > cosines[i*k + cIndx])
                cIndx = j;
        }
        data->assignCluster(i, cIndx);
        // TODO - assignCluster doesn't finalize the assignment, which must
        //        be done in the online case
        // Make another function: assignClusterOnline

        // if the document moved, we must update quality and concept vectors
        // TODO - this part needs to be thread safe (locking)
        if(data->p_asgns[i] != cIndx) {
            // TODO - update the concept vector (live)
            // TODO - all of these updates should only take nz(i) time.
/*            sU[cIndx] = sU[cIndx] - vec_sum(vec_squared(concepts[cIndx]));
            for(int a=0; a<wc; a++) {
                concepts[cIndx][a] =
                    concepts[cIndx][a] + LR * cNorms[cIndx] + doc_matrix[i][a];
            }
            sU[cIndx] = sU[cIndx] + vec_sum(vec_squared(concepts[cIndx]));
            cNorms[cIndx] = sqrt(sU[cIndx]);*/
            
            // remove quality from old vector
            // TODO - ClusterData stores documents in compressed format
            /*float old_quality = vec_dot(data->doc_matrix[i],
                                        data->concepts[data->p_asgns[i]],
                                        data->wc);
            data->qualities[data->p_asgns[i]] -= old_quality;
            // add quality to new vector
            float new_quality = vec_dot(data->doc_matrix[i],
                                        data->concepts[cIndx],
                                        data->wc);
            data->qualities[cIndx] += new_quality;*/
        }
    }
};



// Priority computation for the online version of the algorithm.
struct ComputePriority {
    unsigned int operator() (const int& document_index) const {
        // TODO - implement priority scheme and reference to documents
        return 1;
    }
};



// Run the spherical K-means algorithm using the Galois library.
ClusterData* SPKMeansGalois::runSPKMeans()
{
    /*// first, convert the document matrix to a graph
    Galois::Graph::LC_CSR_Graph<DataNode, float> g;
    
    // set up the Galois structures
    int num_nodes = dc + wc;
    int num_edges = wc;*/


    //-----------------------------------------------------------------------//

    // apply the TXN scheme on the document vectors (normalize them)
    txnScheme();

    // initialize the data arrays; keep track of the arrays locally
    ClusterData *data = new ClusterData(k, dc, wc, doc_matrix);
    float **concepts = data->concepts;
    bool *changed = data->changed;

    // compute initial partitioning, concepts, and quality
    initClusters(data);
    float quality = computeQ(data);
    cout << "Initial quality: " << quality << endl;


    // set up Galois computing structures, and worklist prioritization
    ComputeClustersBasic comp(data);

    // bind the cosineSimilarity function
    comp.cosineSimilarity = bind(&SPKMeans::cosineSimilarity, this,
        placeholders::_1, placeholders::_2, placeholders::_3);

    // this is the worklist ordering scheme using the ComputePriority struct
    // TODO - the ChunkedFIFO: 32 vs. 64 vs. 16 vs. 8 etc.? What does it mean?
    typedef Galois::WorkList::OrderedByIntegerMetric
            <ComputePriority, Galois::WorkList::ChunkedFIFO<32>>
            comp_wl;

    // set up iterators for use by the Galois loops
    auto start_any = boost::make_counting_iterator<int>(0);
    auto end_dc = boost::make_counting_iterator<int>(dc);
    //auto end_k = boost::make_counting_iterator<int>(k);


    // keep track of the run time of this algorithm
    Timer timer;
    timer.start();

    // keep track of all individual component times for analysis
    Timer ptimer;
    Timer ctimer;

    // do spherical k-means loop
    float dQ = Q_THRESHOLD * 10;
    int iterations = 0;
    while(dQ > Q_THRESHOLD) {
        iterations++;

        // compute new partitions based on old concept vectors
        ptimer.start();
        Galois::for_each(start_any, end_dc, comp,
                         Galois::wl<comp_wl>(),
                         Galois::loopname("Compute Clusters"));
        ptimer.stop();

        if(optimize)
            data->findChangedClusters();
        data->applyAssignments();

        // compute new concept vectors and quality
        ctimer.start();
        float n_quality = computeConcepts(data);
        dQ = n_quality - quality;
        quality = n_quality;
        ctimer.stop();

        // report the quality of the current partitioning
        reportQuality(data, quality, dQ);
    }


    // report runtime statistics
    timer.stop();
    reportTime(iterations, timer.get(), ptimer.get(), ctimer.get());

    // return the resulting partitions and concepts in the ClusterData struct
    return data;
}
