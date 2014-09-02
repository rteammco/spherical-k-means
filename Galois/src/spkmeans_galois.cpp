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



// GLOBAL VARIABLES
unsigned int k;
unsigned int dc;
unsigned int wc;
unsigned int nz;



/*********************** GALOIS STRUCTS AND DEFINITIONS ***********************/


// Graph node data defined here (for both documents and features).
struct Node {
    // TODO - fill this out appropriately
    unsigned int node_id;
    unsigned int cluster_id;
    float *cosines;
};


// Define this graph type as just "Graph" for convenience.
typedef Galois::Graph::LC_CSR_Graph<Node, float> Graph;


// The core processing struct used by Galois to run SPKMeans.
struct PartitioningBasic {

    Graph *g;

    float **concepts;
    float *qualities;
    bool *changed;

    // Keep track of the Graph object, and initialize all data arrays.
    PartitioningBasic(Graph &g_)
    {
        g = &g_;

        concepts = new float*[k];
        qualities = new float[k];
        changed = new bool[k];
        for(int i=0; i<k; i++) {
            changed[i] = true;
        }
    }


    // Delete all arrays to clean up memory.
    ~PartitioningBasic()
    {
        delete[] concepts; // TODO - does this delete the 2D list?
        delete[] qualities;
        delete[] changed;
    }


    // Computes the cosine similarity value between the given document and
    // the concept vector identified by the cluster_id.
    // TODO - pass as reference, or no?
    float cosineSimilarity(Graph::GraphNode node, unsigned int cluster_id)
    {
        // sum up the cosine similarity (for each edge)
        float cos = 0;
        for(Graph::edge_iterator itr = g->edge_begin(node, Galois::NONE),
                                 end = g->edge_end(node, Galois::NONE);
                                 itr != end; itr++)
        {
            Graph::GraphNode word = g->getEdgeDst(itr);
            unsigned int word_id = (word - dc);
            float val = g->getEdgeData(itr, Galois::NONE);
            cos += concepts[cluster_id][word_id] * val;
        }

        return cos;
    }


    // Spherical k-means clustering step: assigns the document node to the
    // cluster of closest cosine similarity.
    void operator () (Graph::GraphNode node,
                      Galois::UserContext<Graph::GraphNode> &ctx)
    {
        Node &data = g->getData(node);

        int cIndx = 0;
        if(changed[0])
            data.cosines[0] = cosineSimilarity(node, 0);
        for(int j=0; j<k; j++) {
            if(changed[j])
                data.cosines[j] = cosineSimilarity(node, j);
            if(data.cosines[j] > data.cosines[cIndx])
                cIndx = j;
        }

        unsigned int last_cluster = data.cluster_id;
        data.cluster_id = cIndx;
        if(cIndx != last_cluster) {
            changed[cIndx] = true;
            changed[last_cluster] = true;
        }
    }

};


/******************************************************************************/



void reportQuality(PartitioningBasic &p, float quality, float dQ)
{
    cout << "Quality: " << quality << " (+" << dQ << ")";
    int num_same = 0;
    for (int i=0; i<k; i++)
        if(!(p.changed[i]))
            num_same++;
    cout << " --- " << num_same << " clusters are the same." << endl;
}

float computeQ(PartitioningBasic &p, Graph &g)
{
    float quality = 0;

    //float **sums = new float[k][wc] <- zeros;
    //for each document:
    //  cIndx = document.cluster_id;
    //  for each word in document:
    //      sums[cIndx][word] += word.val;
    //for each i to k:
    //  if changed[i]:
    //      concepts[i] = vec_normalize(sums[i], wc);
    //      qualities[i] = vec_dot(sums[i], concepts[i], wc);
    //TODO - double-check this algorithm... it looks fishy!

    /*for(int i=0; i<k; i++) {
        if(p.changed[i]) {
            float sum_p[wc];
            for(int j=0; j<wc; j++)
                sum_p[j] = 0;
            // add all documents associated with this cluster
            for(int j=0; j<dc; j++) {
                if(data.p_asgns[j] == i) {
                    for(int a=0; a<wc; a++)
                        sum_p[a] += doc_matrix[j][a];
                }
            }
            data->qualities[i] = vec_dot(sum_p, data->concepts[i], wc);
            delete[] sum_p;
        }
        quality += p.qualities[i];
    }*/

    return quality;
}



// Initializes the graph data to prepare it for clustering.
// Returns the number of document nodes in the graph.
void init(Graph &g)
{
    unsigned int max_id = 0;
    unsigned int num_nz = 0;

    for (Graph::iterator it = g.begin(), end = g.end(); it != end; it++) {
        Graph::GraphNode node = *it;
        Node &data = g.getData(node);

        unsigned int num_edges =
            g.edge_end(node, Galois::NONE) - g.edge_begin(node, Galois::NONE);
        num_nz += num_edges;

        // if this node has edges, then it is a document node, so set max_id
        unsigned int cur_id = (unsigned int) node;
        if (num_edges > 0 && max_id < cur_id)
            max_id = cur_id;

        // if document, initialize its node data
        if (num_edges > 0) {
            data.node_id = cur_id;
            data.cluster_id = 0; // TODO
            data.cosines = new float[k];
        }
    }

    // initialize global graph data values
    dc = max_id;
    wc = g.size() - dc;
    nz = num_nz;
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
    k = atoi(argv[2]);
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


    // initialize the graph data and get initial quality
    init(g);
    cout << "Data: " << dc << " documents, " << wc << " words, and "
                     << nz << " non-zero features." << endl;
    PartitioningBasic p(g);
    float quality = 0;//computeQ(p);


    // run spherical k-means
    const float Q_THRESHOLD = 0.001;
    float dQ = Q_THRESHOLD * 10;
    int iterations = 0;
    Galois::Timer timer;
    timer.start();
    while(dQ > Q_THRESHOLD) {
        iterations++;

        // TODO - do_all instead?
        Galois::for_each(g.begin(), g.begin() + dc, p,
                         Galois::loopname("Partitioning"));
                         // Galois::wl...
        
        // TODO - galois version (or just online version)
        for(int i=0; i<k; i++) {
            if(p.changed[i]) {
                //delete[] p.concepts[i];
                //p.concepts[i] = computeConcept();
            }
        }

        // TODO - galois version (or just online version)
        float n_quality = 0;//computeQ(p);
        dQ = n_quality - quality;
        quality = n_quality;
    }


    return 0;
}
