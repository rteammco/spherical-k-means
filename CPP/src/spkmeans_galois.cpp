/* File: spkmeans_galois.cpp
 *
 * Defines the Galois version of the SKPMeans class. This version directly
 * uses the Galois parallel library for multithreaded computation to speed
 * up the algorithm's runtime on multicore systems.
 */

#include "spkmeans.h"

#include <functional>
#include <iostream>

#include "Galois/Galois.h"
#include "Galois/Timer.h"
//#include "Galois/Graph/Graph.h"
//#include "Galois/Graph/LCGraph.h"
#include "llvm/ADT/SmallVector.h"
#include "Galois/LargeArray.h"

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



/****** GALOIS IMPLEMENTATION STARTS HERE ******/


// General compute struct contains variables and the parallel execution method.
/*struct Compute {

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
*/

/*struct TestCompute {
    int wc;
    void operator () (float* &x, Galois::Runtime::UserContextAccess<float*>::SuperTy &ctx) {
        float z = 0;
        for(int i=0; i<wc; i++)
            z += x[i];
        cout << z << endl;
    }
};*/

struct ComputePartitions {

    ClusterData *data;
    float **doc_matrix;

    // just generate the new partitions locally, and swap memory after
    //vector<float*> *new_partitions;
    //Galois::LargeArray<float*> *new_partitions;
    llvm::SmallVector<float*, 5> *new_partitions; // TODO - what the heck is N (5)?

    // binding for the cosineSimilarity function in the SPKMeans object
    function<float(float*, int)> cosineSimilarity;

    // constructor: get the cluster data pointers
    ComputePartitions(ClusterData *data_, float **doc_matrix_)
        : data(data_), doc_matrix(doc_matrix_) {}

    void prepare()
    {
        //new_partitions = new vector<float*>[data->k];
        //new_partitions = new Galois::LargeArray<float*>[data->k];
        new_partitions = new llvm::SmallVector<float*, 5>[data->k];
    }

    // Compute the new partitions:
    void operator () (int &i, Galois::UserContext<int> &ctx)
        //Galois::Runtime::UserContextAccess<float*>::SuperTy &ctx)
    {
        // get pointers to everything from the data struct
        int k = data->k;
        float **concepts = data->concepts;
        bool *changed = data->changed;
        float *cValues = data->cValues;

        // find the best partition for document i
        int cIndx = 0;
        if(changed[0])
            cValues[i*k] = cosineSimilarity(concepts[0], i);
        for(int j=1; j<k; j++) {
            if(changed[j])
                cValues[i*k + j] = cosineSimilarity(concepts[j], i);
            if(cValues[i*k + j] > cValues[i*k + cIndx])
                cIndx = j;
        }

        // add the document to the partition (race condition region)
        new_partitions[cIndx].push_back(doc_matrix[i]);
    }
};


// Run the spherical K-means algorithm using the Galois library.
ClusterData* SPKMeansGalois::runSPKMeans()
{
    // keep track of the run time of this algorithm
    Galois::Timer timer;
    timer.start();

    // apply the TXN scheme on the document vectors (normalize them)
    txnScheme();

    // initialize the data arrays
    ClusterData *data = new ClusterData(k, dc, wc);

    // choose an initial partitioning, and get first concepts
    initPartitions(data);


    // create the Galois computation structs
    ComputePartitions cP(data, doc_matrix);

    // bind the cosine similarity function from SPKMeans to the struct
    cP.cosineSimilarity = bind(&SPKMeans::cosineSimilarity,
                               this, // use this object's instance variables
                               placeholders::_1, placeholders::_2);
    

    // compute initial quality
    float quality = computeQ(data);
    cout << "Initial quality: " << quality << endl;

    // set up iterators for use by the Galois loops
    auto start_dc = boost::make_counting_iterator<int>(0);
    auto end_dc = boost::make_counting_iterator<int>(dc);

    cP.prepare();
    Galois::for_each(start_dc, end_dc, cP, Galois::loopname("Compute Partitions"));
    for(int i=0; i<k; i++)
        cout << cP.new_partitions[i].size() << endl;


    // convert data to Galois format
    /*Galois::LargeArray<float*> docs;
    docs.create(dc, nullptr);

    for(int i=0; i<dc; i++)
        docs[i] = doc_matrix[i];*/

    /*TestCompute t;
    t.wc = wc;
    Galois::for_each(docs.begin(), docs.begin() + dc,
                     t,
                     Galois::loopname("Compute New Partitions"));
    */

    return data;

    // initialize the data arrays
/*    ClusterData *data = new ClusterData(k, dc, wc);

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
        // TODO
        //Galois::for_each(g.begin(), g.begin() + dc,
        //                 cP(),
        //                 Galois::loopname("Compute New Partitions"));
        ptimer.stop();
        p_time += ptimer.get();

        ctimer.start();
        // TODO
        //Galois::for_each(g.begin(), g.begin() + dc,
        //                 cC(),
        //                 Galois::loopname("Compute Concept Vectors"));
        ctimer.stop();
        c_time += ctimer.get();

        qtimer.start();
        // TODO
        //Galois::for_each(g.begin(), g.begin() + dc,
        //                 cQ(),
        //                 Galois::loopname("Compute Quality"));
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
    */
}
