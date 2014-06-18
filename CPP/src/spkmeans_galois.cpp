/* File: spkmeans_galois.cpp
 */

#include "spkmeans.h"

#include <iostream>

#include "Galois/Galois.h"
#include "Galois/Timer.h"
#include "Galois/Graph/Graph.h"
#include "Galois/Graph/LCGraph.h"

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



/****** GALOIS IMPLEMENTATION STARTS HERE ******/

// Define the Galois graph as needed: Nodes are just document vectors.
typedef float* Node;
typedef Galois::Graph::LC_CSR_Graph<Node, void> Graph;



// Galois struct contains variables and the parallel execution method.
struct Compute {

    ClusterData *data;
    Compute(ClusterData *data) : data(data) {}

    void operator () (Graph::GraphNode document,
                      Galois::UserContext<Graph::GraphNode> &ctx) {}

};


struct ComputePartitions : public Compute {

    ComputePartitions(ClusterData *data) : Compute(data) {}

    void operator () (Graph::GraphNode document,
                      Galois::UserContext<Graph::GraphNode> &ctx)
    {
        // compute partitons here
    }
    
};

struct ComputeConcepts : Compute {

    ComputeConcepts(ClusterData *data) : Compute(data) {}

    void operator () (Graph::GraphNode document,
                      Galois::UserContext<Graph::GraphNode> &ctx)
    {
        // compute concept vectors here
    }
    
};

struct ComputeQuality : Compute {

    ComputeQuality(ClusterData *data) : Compute(data) {}

    void operator () (Graph::GraphNode document,
                      Galois::UserContext<Graph::GraphNode> &ctx)
    {
        // compute partition quality here 
    }
    
};

ClusterData* SPKMeansGalois::runSPKMeans()
{
    // keep track of the run time of this algorithm
    Galois::Timer timer;
    timer.start();

    // apply the TXN scheme on the document vectors (normalize them)
    txnScheme();

    // initialize the data arrays
    ClusterData *data = new ClusterData(k, dc, wc);

    // choose an initial partitioning, and get concepts and quality
    initPartitions(data);
    float quality = computeQ(data->partitions, data->p_sizes, data->concepts);

    // initialize the Galois computational components
    ComputePartitions cP(data);
    ComputeConcepts cC(data);
    ComputeQuality cQ(data);

    // TODO - init graph; how?
    Graph g;
    //g.add(???);

    // keep track of all individual component times for analysis
    Galois::Timer ptimer;
    Galois::Timer ctimer;
    Galois::Timer qtimer;
    float p_time = 0;
    float c_time = 0;
    float q_time = 0;

    // do spherical k-means loop
    float dQ = Q_THRESHOLD * 10;
    int iterations = 0;
    while(dQ > Q_THRESHOLD) {
        iterations++;

        // TODO - better to do 3 steps, or all at the same time with Galois?

        ptimer.start();
        /* TODO
        Galois::for_each(g.begin(), g.begin() + dc,
                         cP(),
                         Galois::loopname("Compute New Partitions"));
        */
        ptimer.stop();
        p_time += ptimer.get();

        ctimer.start();
        /* TODO
        Galois::for_each(g.begin(), g.begin() + dc,
                         cC(),
                         Galois::loopname("Compute Concept Vectors"));
        */
        ctimer.stop();
        c_time += ctimer.get();

        qtimer.start();
        /* TODO
        Galois::for_each(g.begin(), g.begin() + dc,
                         cQ(),
                         Galois::loopname("Compute Quality"));
        */
        qtimer.stop();
        q_time += qtimer.get();

        // TODO
        dQ = 0;
        cout << "Quality: " << quality << " (+" << dQ << ")" << endl;
    }

    // report runtime statistics
    timer.stop();
    reportTime(0, timer.get());

    // return the resulting partitions and concepts in the ClusterData struct
    return data;
}
