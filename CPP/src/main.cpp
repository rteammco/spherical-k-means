/* File: spkmeans.cpp
 *
 * A C++ implementation of the Spherical K-means algorithm with an emphasis
 * on parallel computing. This code features three modules: a normal,
 * single-thread implementation, a multithreaded implementation using OpenMP,
 * and a multithreaded implementation using the Galois library:
 *  (http://iss.ices.utexas.edu/?p=projects/galois).
 * The purpose of this code is to compare and optimize the algorithm with
 * different parallel paradigms.
 */


// PROGRAM VERSION
#define VERSION "0.1 (dev)"


#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

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



// Prints a message on how to use this program.
void printUsage()
{
    cout << "Argument options:" << endl
         << "  [-d docfile]     document file path" << endl
         << "  [-v vocabfile]   vocabulary file name" << endl
         << "  [-k num]         value of k (number of clusters)" << endl
         << "  [-t numthreads]  number of threads (when applicable)" << endl
         << "  [--galois]       run in Galois mode" << endl
         << "  [--openmp]       run in OpenMP mode" << endl
         << "  [--noresults]    squelch results from being printed" << endl
         << "  [--autok]        set K automatically using input data" << endl
         << "Other commands:" << endl
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



/* Takes argc and argv from program input and parses the parameters to set
 * values for the k-means algorithm.
 * PARAMETERS (will be filled in be dereferencing except argc and argv):
 *  argc, argv   - The program arguments as passed into the executable.
 *  doc_fname    - String pointer to contain the document file name.
 *  vocab_fname  - String pointer to contain the vocabulary file name.
 *  k            - Int pointer that will be filled with the size of k.
 *  num_threads  - Int pointer that will be filled with the number of threads.
 *  run_type     - Int pointer to indicate which module to run (e.g. OpenMP).
 *  show_results - Bool flag to swith displaying results on or off.
 *  auto_k       - Bool flag to switch choosing K automatically on or off.
 * RETURNS:
 *  RETURN_HELP     to print the program help message and exit.
 *  RETURN_VERSION  to print the program version and exit.
 *  RETURN_ERROR    to print the program help message and exit with an error.
 *  RETURN_SUCCESS  to continue normally.
 */
int processArgs(int argc, char **argv,
    string *doc_fname, string *vocab_fname,
    unsigned int *k, unsigned int *num_threads,
    unsigned int *run_type, bool *show_results, bool *auto_k)
{
    // set defaults before proceeding to check arguments
    *doc_fname = DEFAULT_DOC_FILE;
    *vocab_fname = "";
    *k = DEFAULT_K;
    *num_threads = DEFAULT_THREADS;
    *run_type = RUN_NORMAL;
    *show_results = true;
    *auto_k = false;

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

        // or if K should be selected automatically
        else if(arg == "--autok" || arg == "--auto")
            *auto_k = true;

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
    bool show_results, auto_k;
    int retval = processArgs(argc, argv, &doc_fname, &vocab_fname,
        &k, &num_threads, &run_type, &show_results, &auto_k);
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

    // read data from the document file
    int dc, wc, non_zero;
    float **D = readDocFile(doc_fname.c_str(), &dc, &wc, &non_zero);
    cout << "DATA: " << dc << " documents, " << wc << " words ("
         << non_zero << " non-zero entries)." << endl;

    // if choosing K was set automatically, compute that here if possible
    if(auto_k) {
        int numerator = dc * wc;
        if(non_zero <= 0 || non_zero > numerator)
            cout << "Could not set K automatically. Using k=" << k << endl;
        else
            k = numerator / non_zero;
    }
    cout << "Running SPK Means on \"" << doc_fname << "\" with k=" << k;

    // run the program based on the run type provided (none, openmp, galois)
    ClusterData *data = 0;
    if(run_type == RUN_GALOIS) {
        // tell Galois the max thread count
        SPKMeansGalois spkm(num_threads);
        cout << " [Galois: " << spkm.getNumThreads() << " threads]." << endl;
        data = spkm.runSPKMeans(D, k, dc, wc);
    }
    else if(run_type == RUN_OPENMP) {
        // tell OpenMP the max thread count
        SPKMeansOpenMP spkm(num_threads);
        cout << " [OpenMP: " << num_threads << " threads]." << endl;
        data = spkm.runSPKMeans(D, k, dc, wc);
    }
    else {
        SPKMeans spkm;
        cout << " [single thread]." << endl;
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
