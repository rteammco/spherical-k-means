/* File: spkmeans_galois.cpp
 *
 * Defines the Galois version of the SKPMeans class. This version directly
 * uses the Galois parallel library for multithreaded computation to speed
 * up the algorithm's runtime on multicore systems.
 */

#include "spkmeans.h"

#include <functional>
#include <iostream>
#include <mutex>

#include "Galois/Galois.h"
#include "Galois/Graph/Graph.h"
#include "Galois/Timer.h"
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
        t_ = 2; // TODO - no support in API for this?
    Galois::setActiveThreads(t_);
    num_threads = Galois::getActiveThreads();
}


// Returns the actual number of threads that Galois will use.
unsigned int SPKMeansGalois::getNumThreads()
{
    return num_threads;
}



/*** GALOIS PROCESSING STRUCTURES ***/
struct DocumentNode {
    unsigned int cluster;
};

typedef Galois::Graph::LC_CSR_Graph<DocumentNode, float> Graph;
struct SPKMeansIteration {

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

};



// Run the spherical K-means algorithm using the Galois library.
ClusterData* SPKMeansGalois::runSPKMeans()
{
    // first, convert the document matrix to a graph
    Galois::Graph::LC_CSR_Graph<DocumentNode, float> g;
    

    // keep track of the run time of this algorithm
    Galois::Timer timer;
    timer.start();

    // apply the TXN scheme on the document vectors (normalize them)
    txnScheme();

    // set up the Galois structures
    int num_nodes = dc + wc;
    int num_edges = wc;

    // initialize the data arrays
    ClusterData *data = new ClusterData(k, dc, wc, doc_matrix);

    // choose an initial partitioning, and get first concepts
    initClusters(data);

    /*
    // create the Galois computation structs, and bind the necessary functions
    ComputePartitions cP(data, doc_matrix);
    cP.cosineSimilarity = bind(&SPKMeans::cosineSimilarity,
                               this, // use this object's instance variables
                               placeholders::_1, placeholders::_2);
    ComputeConcepts cC(data);
    cC.computeConcept = bind(&SPKMeans::computeConcept,
                             this, // use this object's instance variables
                             placeholders::_1, placeholders::_2);
    

    // keep track of all individual component times for analysis
    Galois::Timer ptimer;
    Galois::Timer ctimer;
    Galois::Timer qtimer;
    float p_time = 0;
    float c_time = 0;
    float q_time = 0;

    // compute initial quality
    float quality = computeQ(data);
    cout << "Initial quality: " << quality << endl;

    // set up iterators for use by the Galois loops
    auto start_any = boost::make_counting_iterator<int>(0);
    auto end_dc = boost::make_counting_iterator<int>(dc);
    auto end_k = boost::make_counting_iterator<int>(k);

    // do spherical k-means loop
    float dQ = Q_THRESHOLD * 10;
    int iterations = 0;
    while(dQ > Q_THRESHOLD) {
        iterations++;

        // compute new partitions based on old concept vectors
        ptimer.start();
        cP.prepare();
        Galois::for_each(start_any, end_dc, cP,
            Galois::loopname("Compute Partitions"));
        ptimer.stop();
        p_time += ptimer.get();

        // check if partitions changed since last time
        findChangedPartitionsUnordered(cP.new_partitions, data);

        // transfer new partitions to the partitions array
        copyPartitions(cP.new_partitions, data);

        // compute new concept vectors
        ctimer.start();
        Galois::for_each(start_any, end_k, cC,
            Galois::loopname("Compute Concepts"));
        ctimer.stop();
        c_time += ctimer.get();

        // compute quality of new partitioning
        qtimer.start();
        float n_quality = computeQ(data);
        dQ = n_quality - quality;
        quality = n_quality;
        qtimer.stop();
        q_time += qtimer.get();

        // report quality and (if optimizing) which partitions changed
        reportQuality(data, quality, dQ);
    }


    // report runtime statistics
    timer.stop();
    reportTime(iterations, timer.get(), p_time, c_time, q_time);

    // clean up
    cP.clearMemory();
    */

    // return the resulting partitions and concepts in the ClusterData struct
    return data;
}
