/* File: spkmeans_galois.cpp
 *
 * A pure Galois-only implementation of spherical k-means. This is a simpler
 * non-object-oriented version that reads data directly from Galois graph
 * files. It is a "prototype" implementation.
 */

#include <iostream>
#include <string>

#include "Galois/Galois.h"
#include "Galois/Graph/Graph.h"
//#include "Galois/LCGraph.h"
#include "Galois/Timer.h"
//#include "Galois/Statistic.h"

#include "boost/iterator/counting_iterator.hpp"

using namespace std;



/*********************** GALOIS STRUCTS AND DEFINITIONS ***********************/


// Graph node data defined here (for both documents and features).
struct Node {
    // TODO - fill this out appropriately
    unsigned int cluster_id;
};


// Define this graph type as just "Graph" for convenience.
typedef Galois::Graph::LC_CSR_Graph<Node, float> Graph;


// The core processing struct used by Galois to run SPKMeans.
struct PartitionBasic {
    ;
};


/******************************************************************************/



// Initializes the graph data to prepare it for clustering.
// Returns the number of document nodes in the graph.
unsigned int init(Graph &g)
{
    unsigned int max_id = 0;

    for (Graph::iterator it = g.begin(), end = g.end(); it != end; it++) {
        Graph::GraphNode node = *it;
        Node &data = g.getData(node);

        unsigned int num_edges =
            g.edge_end(node, Galois::NONE) - g.edge_begin(node, Galois::NONE);

        // if this node has edges, then it is a document node, so set max_id
        unsigned int cur_id = (unsigned int) node;
        if (num_edges > 0 && max_id < cur_id)
            max_id = cur_id;

        // if document, initialize its node data
        if (num_edges > 0) {
            // TODO - fill this out appropriately
            data.cluster_id = 0;
        }
    }

    return max_id;
}



// Main (start here):
// Reads arguments and graph file, and sets up Galois.
int main(int argc, char ** argv)
{
    // process arguments parameters
    if (argc < 3) {
        cout << "Usage: <.gr data file> <k> <threads>" << endl;
        return 0;
    }
    const char *file = argv[1];
    const unsigned int k = atoi(argv[2]);
    unsigned int n_threads = atoi(argv[3]);

    // initialize Galois threads
    Galois::setActiveThreads(n_threads);
    n_threads = Galois::getActiveThreads();
    string s = "s";
    if (n_threads == 1)
        s = "";
    cout << "Running with " << n_threads << " thread" << s << "." << endl;

    // read graph file
    Graph g;
    Galois::Graph::readGraph(g, file);

    // initialize the graph data
    unsigned int dc = init(g);
    unsigned int wc = g.size() - dc;
    cout << "Data: " << dc << " documents and " << wc << " words." << endl;

    return 0;
}
