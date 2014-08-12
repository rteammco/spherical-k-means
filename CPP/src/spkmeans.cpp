/* File: spkmeans.cpp
 *
 * Defines the basic SPKMeans functions for the generic SPKMeans class.
 * This class provides a single-threaded implementation of the algorithm,
 * as well as some basic functions useful for the multithreaded versions.
 */

#include "spkmeans.h"

#include "timer.h"
#include "vectors.h"

#include <iostream>
#include <vector>

using namespace std;



// Constructor: initialize variables and document norms.
SPKMeans::SPKMeans(float **doc_matrix_, int k_, int dc_, int wc_)
    : doc_matrix(doc_matrix_), k(k_), dc(dc_), wc(wc_)
{
    // init the doc vector norms
    doc_norms = new float[dc];
    for(int i=0; i<dc; i++)
        doc_norms[i] = vec_norm(doc_matrix[i], wc);

    // default scheme to TXN
    prep_scheme = TXN_SCHEME;

    // initialize optimization to true (optimization will happen)
    optimize = true;
}



// Destructor: clean up document norms array.
SPKMeans::~SPKMeans()
{
    delete[] doc_norms;
}



// Set the current scheme to the given type.
void SPKMeans::setScheme(SPKMeans::Scheme type)
{
    prep_scheme = type;
}



// Disables optimization for testing purposes.
void SPKMeans::disableOptimization()
{
    optimize = false;
}



// Enables optimization if it was disabled.
void SPKMeans::enableOptimization()
{
    optimize = true;
}



// Reports the overall quality and, if optimizing, also displays how many
// clusters have changed.
void SPKMeans::reportQuality(ClusterData *data, float quality, float dQ)
{
    cout << "Quality: " << quality << " (+" << dQ << ")";
    if(optimize) {
        int num_same = 0;
        for (int i=0; i<k; i++)
            if(!(data->changed[i]))
                num_same++;
        cout << " --- " << num_same << " clusters are the same." << endl;
    }
    else
        cout << " --- optimization disabled." << endl;
}



// Reports time data after running the algorithm.
void SPKMeans::reportTime(int iterations, float total_time,
                          float p_time, float c_time, float q_time)
{
    cout << "Done in " << total_time / 1000
         << " seconds after " << iterations << " iterations." << endl;
    float total = p_time + c_time + q_time;
    if(total == 0)
        cout << "No individual time stats available." << endl;
    else {
        cout << "Timers (ms): " << endl
             << "   partitioning [" << p_time << "] ("
                << (p_time/total)*100 << "%)" << endl
             << "   concepts     [" << c_time << "] ("
                << (c_time/total)*100 << "%)" << endl
             << "   quality      [" << q_time << "] ("
                << (q_time/total)*100 << "%)" << endl;
    }
}



// Applies the TXN scheme to each document vector of the given matrix.
// TXN effectively just normalizes each of the document vectors.
// TXN scheme must be set, otherwise this function will do nothing.
void SPKMeans::txnScheme()
{
    if(prep_scheme != SPKMeans::TXN_SCHEME)
        return;

    for(int i=0; i<dc; i++)
        vec_normalize(doc_matrix[i], wc);
}



// Initializes the first partitioning (randomly or otherwise assigned) to
// provide a starting point for the clustering algorithm.
void SPKMeans::initClusters(ClusterData *data)
{
    // choose an initial partitioning
    int split = dc / k;
    cout << "Split = " << split << endl;
    int base = 1;
    for(int i=0; i<k; i++) {
        int top = base + split - 1;
        if(i == k-1)
            top = dc;
        int p_size = top - base + 1;
        for(int j=0; j<p_size; j++)
            data->p_asgns[base-1 + j] = i;
        base = base + split;
    }

    // compute the initial concept vectors
    for(int i=0; i<k; i++)
        data->concepts[i] = computeConcept(data, i);
}



// Returns the total quality of all clusters by summing the qualities of
// each individual cluster. If optimization is enabled, uses cached values
// whenever possible.
float SPKMeans::computeQ(ClusterData *data)
{
    float quality = 0;

    // NOTE - there was a single-cluster computeQ function, re-implement?
    for(int i=0; i<k; i++) {
        if(data->changed[i]) {
            float *sum_p = vec_zeros(wc);
            // add all documents associated with this cluster
            for(int j=0; j<dc; j++) {
                if(data->p_asgns[j] == i) {
                    for(int a=0; a<wc; a++)
                        sum_p[a] += doc_matrix[j][a];
                }
            }
            data->qualities[i] = vec_dot(sum_p, data->concepts[i], wc);
            delete[] sum_p;
        }
        quality += data->qualities[i];
    }

    return quality;
}



// Computes the cosine similarity value of the two given vectors (dv and cv).
float SPKMeans::cosineSimilarity(ClusterData *data, int doc_index, int cIndx)
{
    // TODO - same optimization for concepts? can we cache doc norms better?
    float cnorm = vec_norm(data->concepts[cIndx], wc);
    float dnorm = doc_norms[doc_index];

    // here is where we save time: compute the dot product!
    float dotp = 0;
    for(int i=0; i<(data->docs[doc_index].count); i++) {
        int word = data->docs[doc_index].words[i].index;
        float value = data->docs[doc_index].words[i].value;
        dotp += data->concepts[cIndx][word] * value;
    }
    
    return dotp / (dnorm * cnorm);
}



// Computes the concept vector of the given cluster (by index). The
// cluster documents are accessed using the ClusterData struct, and the
// associated concept vector will be allocated and populated.
float* SPKMeans::computeConcept(ClusterData *data, int cIndx)
{
    // create the concept vector and initialize it to 0
    float *concept = new float[wc];
    for(int i=0; i<wc; i++)
        concept[i] = 0;

    // find documents associated w/ this cluster, and sum them
    for(int i=0; i<dc; i++) { // for each document
        if(data->p_asgns[i] == cIndx) { // if doc in cluster
            for(int j=0; j<wc; j++) // add words to concept
                concept[j] += doc_matrix[i][j];
        }
    }

    // normalize the concept vector and return it
    //vec_divide(concept, wc, wc);
    vec_normalize(concept, wc);
    return concept;
}



// Runs the spherical k-means algorithm on the given sparse matrix D and
// clusters the data into k clusters. Non-parallel (standard) version.
ClusterData* SPKMeans::runSPKMeans()
{
    // keep track of the run time for this algorithm
    Timer timer;
    timer.start();

    // keep track of all individual component times for analysis
    Timer ptimer;
    Timer ctimer;
    Timer qtimer;

    // apply the TXN scheme on the document vectors (normalize them)
    txnScheme();

    // initialize the data arrays; keep track of the arrays locally
    ClusterData *data = new ClusterData(k, dc, wc, doc_matrix);
    float **concepts = data->concepts;
    bool *changed = data->changed;
    float *cosines = data->cosine_similarities;

    // compute initial partitioning, concepts, and quality
    initClusters(data);
    float quality = computeQ(data);
    cout << "Initial quality: " << quality << endl;
    

    // do spherical k-means loop
    float dQ = Q_THRESHOLD * 10;
    int iterations = 0;
    while(dQ > Q_THRESHOLD) {
        iterations++;

        // TODO - temporary (testing) -----------------------------------------
        /*if(iterations > 1) {
            cout << "New version." << endl;
            ptimer.start();
            float averagePriority = data->getAverageMovedPriority();
            for(int i=0; i<dc; i++) {
                if(data->doc_priorities[i] > averagePriority) {
                    int cIndx = 0;
                    if(changed[0])
                        cosines[i*k] = cosineSimilarity(data, i, 0);
                    for(int j=1; j<k; j++) {
                        if(changed[j])
                            cosines[i*k + j] = cosineSimilarity(data, i, j);
                        if(cosines[i*k + j] > cosines[i*k + cIndx])
                            cIndx = j;
                    }
                    // compute the priority heuristic and assign the document
                    float priority = 1 - cosines[i*k + data->p_asgns[i]];
                    data->assignCluster(i, cIndx, priority);
                }
            }
        } else {*/
        // --------------------------------------------------------------------

        // compute new clusters based on old concept vectors
        ptimer.start();

        // TODO - temporary testing for empty clusters
        bool has_docs[k];
        for(int i=0; i<k; i++)
            has_docs[i] = false;

        for(int i=0; i<dc; i++) {
            // only update cosine similarities if cluster has changed
            int cIndx = 0;
            if(changed[0])
                cosines[i*k] = cosineSimilarity(data, i, 0);
            for(int j=1; j<k; j++) {
                if(changed[j])
                    cosines[i*k + j] = cosineSimilarity(data, i, j);
                if(cosines[i*k + j] > cosines[i*k + cIndx])
                    cIndx = j;
            }
            // compute the priority heuristic and assign the document
            //float priority = 1 - cosines[i*k + data->p_asgns[i]];
            data->assignCluster(i, cIndx);//, priority);
            has_docs[cIndx] = true;
        }//}

        // TODO - temporary testing for empty clusters:
        for(int i=0; i<k; i++) {
            if(!has_docs[i])
                cout << "Cluster " << i << " is empty!" << endl;
        }

        ptimer.stop();

        // TODO - temporary (testing) -----------------------------------------
        // report priorities and number of documents that moved
        /*cout << "Number of documents moved: " << data->num_moved << endl;
        cout << "Average document priority: "
             << data->getAveragePriority() << endl;
        cout << "   Average moved priority: "
             << data->getAverageMovedPriority() << endl;
        cout << "  Average stayed priority: "
             << data->getAverageStayPriority() << endl;*/
        // --------------------------------------------------------------------

        // update which clusters changed since last time, then swap pointers
        if(optimize)
            data->findChangedClusters();
        data->applyAssignments();

        // compute new concept vectors
        ctimer.start();
        for(int i=0; i<k; i++) {
            if(changed[i]) {
                delete[] concepts[i];
                concepts[i] = computeConcept(data, i);
            }
        }
        ctimer.stop();

        // compute quality of new partitioning
        qtimer.start();
        float n_quality = computeQ(data);
        dQ = n_quality - quality;
        quality = n_quality;
        qtimer.stop();

        // report the quality of the current partitioning
        reportQuality(data, quality, dQ);
    }


    // report runtime statistics
    timer.stop();
    reportTime(iterations, timer.get(), ptimer.get(), ctimer.get(),
               qtimer.get());

    // return the resulting clusters and concepts in the ClusterData struct
    return data;
}
