/* File: spkmeans.cpp
 *
 * A parallel implementation of the Spherical K-Means algorithm using the
 * Galois library (http://iss.ices.utexas.edu/?p=projects/galois).
 */


#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

#include "Galois/Galois.h"
#include "Galois/Timer.h"

#include "cluster_data.h"
#include "reader.h"
#include "vectors.h"


// default parameters
#define DEFAULT_K 2
#define DEFAULT_THREADS 2
#define Q_THRESHOLD 0.001
#define DEFAULT_DOC_FILE "data"


using namespace std;



// Debug: prints the given vector (array) to std out.
void printVec(float *vec, int size)
{
    for(int i=0; i<size; i++)
        cout << vec[i] << " ";
    cout << endl;
}



// Prints a short message on how to use this program.
void printUsage()
{
    cout << "Usage: "
         << "./spkmeans doc_file [k] [num_threads] [vocab_file]"
         << endl;
    cout << "Default values are: k = 2, 2 threads, no vocabulary file."
         << endl;
}



// Applies the TXN scheme to each document vector of the given matrix.
// TXN effectively just normalizes each of the document vectors.
void txnScheme(float **doc_matrix, int dc, int wc)
{
    for(int i=0; i<dc; i++)
        vec_normalize(doc_matrix[i], wc);
}



// Returns the quality of the given partition by doing a dot product against
// its given concept vector.
float computeQuality(float **partition, int p_size, float *concept, int wc)
{
    float *sum_p = vec_sum(partition, wc, p_size);
    float quality = vec_dot(sum_p, concept, wc);
    delete sum_p;
    return quality;
}



// Returns the total quality of all partitions by summing the qualities of
// each individual partition.
float computeQuality(float ***partitions, int *p_sizes, float **concepts,
    int k, int wc)
{
    float quality = 0;
    for(int i=0; i<k; i++)
        quality += computeQuality(partitions[i], p_sizes[i], concepts[i], wc);
    return quality;
}



// Computes the cosine similarity value of the two given vectors (dv and cv).
float cosineSimilarity(float *dv, float *cv, int wc)
{
    return vec_dot(dv, cv, wc) / (vec_norm(dv, wc) * vec_norm(cv, wc));
}



// Computes the concept vector of the given partition. A partition is an array
// of document vectors, and the concept vector will be allocated and populated.
float* computeConcept(float **partition, int p_size, int wc)
{
    float *cv = vec_sum(partition, wc, p_size);
    vec_multiply(cv, wc, (1.0 / wc));
    vec_divide(cv, wc, vec_norm(cv, wc));
    return cv;
}



// Runs the spherical k-means algorithm on the given sparse matrix D and
// clusters the data into k partitions.
ClusterData runSPKMeans(float **doc_matrix, unsigned int k, int dc, int wc)
{
    // keep track of the run time for this algorithm
    Galois::Timer timer;
    timer.start();


    // apply the TXN scheme on the document vectors (normalize them)
    txnScheme(doc_matrix, dc, wc);


    // initialize the data arrays; keep track of the arrays locally
    ClusterData data(k, dc, wc);
    float ***partitions = data.partitions;
    int *p_sizes = data.p_sizes;
    float **concepts = data.concepts;


    // create the first arbitrary partitioning
    int split = floor(dc / k);
    cout << "Split = " << split << endl;
    int base = 1;
    for(int i=0; i<k; i++) {
        int top = base + split - 1;
        if(i == k-1)
            top = dc;

        int p_size = top - base + 1;
        p_sizes[i] = p_size;
        cout << "Created new partition of size " << p_size << endl;

        partitions[i] = new float*[p_size];
        for(int j=0; j<p_size; j++)
            partitions[i][j] = doc_matrix[base + j - 1];

        base = base + split;
    }


    // compute concept vectors
    for(int i=0; i<k; i++)
        concepts[i] = computeConcept(partitions[i], p_sizes[i], wc);


    // compute initial quality of the partitions
    float quality = computeQuality(partitions, p_sizes, concepts, k, wc);
    cout << "Initial quality: " << quality << endl;


    // keep track of all individual component times for analysis
    Galois::Timer ptimer;
    Galois::Timer ttimer;
    Galois::Timer ctimer;
    Galois::Timer qtimer;
    float p_time = 0;
    float t_time = 0;
    float c_time = 0;
    float q_time = 0;

    // do spherical k-means loop
    float dQ = Q_THRESHOLD * 10;
    int iterations = 0;
    while(dQ > Q_THRESHOLD) {
        iterations++;

        // compute new partitions based on old concept vectors
        ptimer.start();
        vector<float*> *new_partitions = new vector<float*>[k];
        for(int i=0; i<k; i++)
            new_partitions[i] = vector<float*>();
        for(int i=0; i<dc; i++) {
            int cIndx = 0;
            float cVal = cosineSimilarity(doc_matrix[i], concepts[0], wc);
            for(int j=1; j<k; j++) {
                float new_cVal = cosineSimilarity(doc_matrix[i], concepts[j], wc);
                if(new_cVal > cVal) {
                    cVal = new_cVal;
                    cIndx = j;
                }
            }
            new_partitions[cIndx].push_back(doc_matrix[i]);
        }
        ptimer.stop();
        p_time += ptimer.get();

        // transfer the new partitions to the partitions array
        ttimer.start();
        data.clearPartitions();
        for(int i=0; i<k; i++) {
            partitions[i] = new_partitions[i].data();
            p_sizes[i] = new_partitions[i].size();
        }
        ttimer.stop();
        t_time += ttimer.get();

        // compute new concept vectors
        ctimer.start();
        data.clearConcepts();
        for(int i=0; i<k; i++)
            concepts[i] = computeConcept(partitions[i], p_sizes[i], wc);
        ctimer.stop();
        c_time += ctimer.get();

        // compute quality of new partitioning
        qtimer.start();
        float n_quality = computeQuality(partitions, p_sizes, concepts, k, wc);
        dQ = n_quality - quality;
        quality = n_quality;
        qtimer.stop();
        q_time += qtimer.get();

        cout << "Quality: " << quality << " (+" << dQ << ")" << endl;
    }


    // report runtime statistics
    timer.stop();
    float time_in_ms = timer.get();
    cout << "Done in " << time_in_ms / 1000
         << " seconds after " << iterations << " iterations." << endl;
    float total = p_time + t_time + c_time + q_time;
    cout << "Timers (ms): " << endl
         << "   partition [" << p_time << "] (" << (p_time/total)*100 << "%)" << endl
         << "   p_transfer [" << t_time << "] (" << (t_time/total)*100 << "%)" << endl
         << "   concepts [" << c_time << "] (" << (c_time/total)*100 << "%)" << endl
         << "   quality [" << q_time << "] (" << (q_time/total)*100 << "%)" << endl;


    // return the resulting partitions and concepts in the ClusterData struct
    return data;
}



// Displays the results of each partition. If a words list is provided,
// the top num_to_show words will be displayed for each partition.
// Otherwise, if words is a null pointer, only the indices will be shown.
void displayResults(ClusterData *data, char **words, int num_to_show = 10)
{
    // make sure num_to_show doesn't exceed the actual word count
    if(num_to_show > data->wc)
        num_to_show = data->wc;

    // for each partition, sum the weights of each word, and show the top
    //  words that occur in the partition:
    for(int i=0; i<(data->k); i++) {
        cout << "Partition #" << (i+1) << ":" << endl;
        // sum the weights
        float *sum = vec_sum(data->partitions[i], data->wc, data->p_sizes[i]);

        // sort this sum using C++ priority queue (keeping track of indices)
        vector<float> values(sum, sum + data->wc);
        priority_queue<pair<float, int>> q;
        for(int i=0; i<values.size(); i++)
            q.push(pair<float, int>(values[i], i));

        // show top num_to_show words
        for(int i=0; i<num_to_show; i++) {
            int index = q.top().second;
            if(words != 0)
                cout << "   " << words[index] << endl;
            else
                cout << "   " << index << endl;
            q.pop();
        }
    }
}



// Takes argc and argv from program input and parses the parameters to set
// values for k (number of clusters) and num_threads (the maximum number
// of threads for Galois to use).
// Returns -1 on fail (provided file doesn't exist), else 0 on success.
int processArgs(int argc, char **argv,
    string *doc_fname, unsigned int *k, unsigned int *num_threads,
    string *vocab_fname)
{
    // set up document file name
    if(argc >= 2)
        *doc_fname = string(argv[1]);
    else
        *doc_fname = DEFAULT_DOC_FILE;

    // check that the file exists - if not, error
    ifstream test(doc_fname->c_str());
    if(!test.good()) {
        cout << "Error: file \"" << *doc_fname << "\" does not exist." << endl;
        test.close();
        return -1;
    }
    test.close();

    // set up size of k
    if(argc >= 3)
        *k = atoi(argv[2]);
    else
        *k = DEFAULT_K;

    // set up number of threads
    if(argc >= 4) 
        *num_threads = atoi(argv[3]);
    else
        *num_threads = DEFAULT_THREADS;

    // set up vocabulary file name
    if(argc >= 5)
        *vocab_fname = argv[4];
    else
        *vocab_fname = "";

    return 0;
}



// main: set up Galois and start the clustering process.
int main(int argc, char **argv)
{
    // get file names, and set up k and number of threads
    string doc_fname, vocab_fname;
    unsigned int k, num_threads;
    if(processArgs(argc, argv,
        &doc_fname, &k, &num_threads, &vocab_fname) != 0)
    {
        printUsage();
        return -1;
    }

    // tell Galois the max thread count
    Galois::setActiveThreads(num_threads);
    num_threads = Galois::getActiveThreads();
    cout << "Running SPK Means on \"" << doc_fname << "\" with k=" << k
         << " (" << num_threads << " threads)." << endl;

    // set up the sparse matrix
    int dc, wc;
    float **D = readDocFile(doc_fname.c_str(), dc, wc);
    cout << dc << " documents, " << wc << " words." << endl;

    // run spherical k-means on the given sparse matrix D
    ClusterData data = runSPKMeans(D, k, dc, wc);

    char **words = readWordsFile(vocab_fname.c_str(), wc);
    displayResults(&data, words, 10);

    return 0;
}
