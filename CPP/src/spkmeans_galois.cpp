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
    // apply txn stream
    // initial clusters

    ClusterData *data = new ClusterData(k, dc, wc);

    ComputePartitions cP(data);
    ComputeConcepts cC(data);
    ComputeQuality cQ(data);

    // init graph
    Graph g;
    //g.add(???);

    return data;
}
