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


// graph node defined here (for both documents and features)
struct Node {
    float *x;
};


// define this graph type as just "Graph"
typedef Galois::Graph::LC_CSR_Graph<Node, float> Graph;


/******************************************************************************/



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
    const int k = atoi(argv[2]);
    int n_threads = atoi(argv[3]);

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

    cout << "It works!" << endl;
    return 0;
}
