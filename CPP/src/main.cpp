/* File: spkmeans.cpp
 *
 * A parallel implementation of the Spherical K-Means algorithm using the
 * Galois library (http://iss.ices.utexas.edu/?p=projects/galois) and OpenMP.
 */

// PROGRAM VERSION
#define VERSION "0.1 (dev)"


#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "Galois/Galois.h"
#include "Galois/Timer.h"

#include "cluster_data.h"
#include "reader.h"
#include "spkmeans.h"
#include "vectors.h"


// return flags (for argument parsing)
#define RETURN_SUCCESS 0
#define RETURN_HELP 1
#define RETURN_VERSION 2
#define RETURN_ERROR -1

// default parameters
#define DEFAULT_K 2
#define DEFAULT_THREADS 2
#define DEFAULT_DOC_FILE "docs"

// type of parallel implementations
#define RUN_NORMAL 0
#define RUN_GALOIS 1
#define RUN_OPENMP 2


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
    // $ ./spkmeans -d docfile -w wordfile -k 2 -t 2 --galois
    cout << "Usage: " << endl
         << " $ ./spkmeans [-d docfile] [-v vocabfile] [-k k] [-t numthreads] "
         << "[--galois OR --openmp]" << endl
         << "Other commands: " << endl
         << " $ ./spkmeans --help" << endl
         << " $ ./spkmeans --version" << endl;
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
    string *doc_fname, string *vocab_fname, unsigned int *k,
    unsigned int *num_threads, unsigned int *run_type, bool *show_results)
{
    // set defaults before proceeding to check arguments
    *doc_fname = DEFAULT_DOC_FILE;
    *vocab_fname = "";
    *k = DEFAULT_K;
    *num_threads = DEFAULT_THREADS;
    *run_type = RUN_NORMAL;
    *show_results = true;

    // check arguments: expected command as follows:
    // $ ./spkmeans -d docfile -w wordfile -k 2 -t 2 --galois
    for(int i=1; i<argc; i++) {
        string arg(argv[i]);

        // if flag is --help, return 1 to print usage instructions
        if(arg == "--help" || arg == "-h")
            return RETURN_HELP;
        // if flag is --version, return 2 to print version number
        else if(arg == "--version" || arg == "-V")
            return RETURN_VERSION;

        // if the flag was to run as galois or openmp, set the run type
        else if(arg == "--galois")
            *run_type = RUN_GALOIS;
        else if(arg == "--openmp")
            *run_type = RUN_OPENMP;

        // also check if flag was set to squelch displaying results
        else if(arg == "--noresults")
            *show_results = false;

        // otherwise, check the given flag value
        else {
            i++;
            if(arg == "-d") // document file
                *doc_fname = string(argv[i]);
            else if(arg == "-w" || arg == "-v") // words file
                *vocab_fname = string(argv[i]);
            else if(arg == "-k") // size of k
                *k = atoi(argv[i]);
            else if(arg == "-t") // number of threads
                *num_threads = atoi(argv[i]);
        }
    }

    // check that the document file exists - if not, return error
    ifstream test(doc_fname->c_str());
    if(!test.good()) {
        cout << "Error: file \"" << *doc_fname << "\" does not exist." << endl;
        test.close();
        return RETURN_ERROR;
    }
    test.close();

    return RETURN_SUCCESS;
}



// main: set up Galois and start the clustering process.
int main(int argc, char **argv)
{
    // get file names, and set up k and number of threads
    string doc_fname, vocab_fname;
    unsigned int k, num_threads, run_type;
    bool show_results;
    int retval = processArgs(argc, argv,
        &doc_fname, &vocab_fname, &k, &num_threads, &run_type, &show_results);
    if(retval == RETURN_ERROR) {
        printUsage();
        return -1;
    }
    else if(retval == RETURN_HELP) {
        printUsage();
        return 0;
    }
    else if(retval == RETURN_VERSION) {
        cout << "Version: " << VERSION << endl;
        return 0;
    }

    // set up the sparse matrix
    int dc, wc;
    float **D = readDocFile(doc_fname.c_str(), dc, wc);
    cout << dc << " documents, " << wc << " words." << endl;
    cout << "Running SPK Means on \"" << doc_fname << "\" with k=" << k;

    // run the program based on the run type provided (none, openmp, galois)
    ClusterData *data = 0;
    /*else if(run_type == RUN_GALOIS) {
        // tell Galois the max thread count
        cout << " [Galois: " << num_threads << " threads]." << endl;
        cout << " (" << num_threads << " threads)." << endl;
        Galois::setActiveThreads(num_threads);
        num_threads = Galois::getActiveThreads();
        data = runSPKMeans(D, k, dc, wc);
    }*/
    if(run_type == RUN_OPENMP) {
        // tell OpenMP the max thread count
        cout << " [OpenMP: " << num_threads << " threads]." << endl;
        SPKMeansOpenMP spkm(num_threads);
        data = spkm.runSPKMeans(D, k, dc, wc);
    }
    else {
        cout << " [single thread]." << endl;
        SPKMeans spkm;
        data = spkm.runSPKMeans(D, k, dc, wc);
    }

    // display the results of the algorithm (if anything happened)
    if(data) {
        if(show_results) {
            char **words = readWordsFile(vocab_fname.c_str(), wc);
            displayResults(data, words, 10);
        }
        delete data;
    }

    return 0;
}
