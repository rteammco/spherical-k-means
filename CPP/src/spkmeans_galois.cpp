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



/****** GALOIS IMPLEMENTATION STARTS HERE ******/

struct ComputePartitions {

    ClusterData *data;
    float **doc_matrix;
    mutex *mut;

    // just generate the new partitions locally, and swap memory after
    //vector<float*> *new_partitions;
    //Galois::LargeArray<float*> *new_partitions;
    llvm::SmallVector<float*, 5> *new_partitions; // TODO - what the heck is N (5)?

    // binding for the cosineSimilarity function in the SPKMeans object
    function<float(float*, int)> cosineSimilarity;

    // constructor: get the cluster data pointers
    ComputePartitions(ClusterData *data_, float **doc_matrix_)
        : data(data_), doc_matrix(doc_matrix_)
    {
        new_partitions = 0;
        mut = new mutex();
    }
    
    // clear memory: delete the mutex
    // can't put in destructor because Galois keeps calling it...?
    void clearMemory()
    {
        delete mut;
    }

    
    // clear out the new_partitions vector for new computations
    // TODO - memory leak: delete the old one
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

        // add the document to the partition (race condition region):
        // - race condition is handled by Galois?
        mut->lock();
        new_partitions[cIndx].push_back(doc_matrix[i]);
        mut->unlock();
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
    auto start_dc = boost::make_counting_iterator<int>(0);
    auto end_dc = boost::make_counting_iterator<int>(dc);

    // do spherical k-means loop
    float dQ = Q_THRESHOLD * 10;
    int iterations = 0;
    while(dQ > Q_THRESHOLD) {
        iterations++;

        // compute new partitions based on old concept vectors
        ptimer.start();
        cP.prepare();
        cout << "Prepare" << endl;
        // TODO - libc error in for_each part
        Galois::for_each(start_dc, end_dc, cP, Galois::loopname("Compute Partitions"));
        ptimer.stop();
        p_time += ptimer.get();
        cout << "Partitions" << endl;

        // check if partitions changed since last time
        //findChangedPartitionsUnordered(cP.new_partitions, data);
        // ^ TODO ^ function works with std::vector, not llvm::SmallVector
        if(optimize) {
            for(int i=0; i<k; i++) {
                if(data->p_sizes[i] == cP.new_partitions[i].size()) {
                    data->changed[i] = false;
                    for(int j=0; j<(data->p_sizes[i]); j++) {
                        if(find(cP.new_partitions[i].begin(),
                                cP.new_partitions[i].end(),
                                data->partitions[i][j]) == cP.new_partitions[i].end())
                        {
                            data->changed[i] = true;
                            break;
                        }
                    }
                }
                else
                    data->changed[i] = true;
            }
        }
        cout << "Change" << endl;

        // transfer new partitions to the partitions array
        //copyPartitions(cP.new_partitions, data);
        // ^ TODO ^ function works with std::vector, not llvm::SmallVector
        data->clearPartitions();
        for(int i=0; i<k; i++) {
            data->partitions[i] = cP.new_partitions[i].data();
            data->p_sizes[i] = cP.new_partitions[i].size();
        }
        cout << "Copy" << endl;

        // compute new concept vectors
        ctimer.start();
        ctimer.stop();
        c_time += ctimer.get();
        cout << "Concepts" << endl;

        // compute quality of new partitioning
        // TODO - segfault somewhere below
        qtimer.start();
        float n_quality = computeQ(data);
        dQ = n_quality - quality;
        quality = n_quality;
        qtimer.stop();
        q_time += qtimer.get();
        cout << "Quality" << endl;

        // report quality and (if optimizing) which partitions changed
        reportQuality(data, quality, dQ);
        cout << "Report" << endl;

        // TODO - temporary
        break;
    }


    // report runtime statistics
    timer.stop();
    reportTime(iterations, timer.get(), p_time, c_time, q_time);

    // clean up
    cP.clearMemory();

    // return the resulting partitions and concepts in the ClusterData struct
    return data;



    // TODO - delete all below

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
