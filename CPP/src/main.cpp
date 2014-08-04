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
#define VERSION "0.2 (dev)"


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
#define DEFAULT_THREADS 0 // 0 means default to max
#define DEFAULT_DOC_FILE "test.txt"

// type of parallel implementations
#define RUN_NORMAL 0
#define RUN_GALOIS 1
#define RUN_OPENMP 2


using namespace std;



// Prints a message on how to use this program.
void printUsage()
{
    cout << endl;
    cout << "$ ./spkmeans [options]" << endl;
    cout << "Argument options:" << endl
         << "  [-d docfile]     set document file path" << endl
         << "  [-v vocabfile]   set vocabulary file path" << endl
         << "  [-k num]         set value of k (number of clusters)" << endl
         << "  [-t numthreads]  set number of threads* (if applicable)" << endl
         << "  [--galois]       run in Galois mode (if available)" << endl
         << "  [--openmp]       run in OpenMP mode" << endl
         << "  [--autok]        set K automatically using input data" << endl
         << "  [--noscheme]     do not normalize weight values" << endl
         << "  [--noresults]    squelch results from being printed" << endl
         << "  [--noop]         turn off all optimizations" << endl
         << "Other commands:" << endl
         << "  $ ./spkmeans --help" << endl
         << "  $ ./spkmeans --version" << endl;
    cout << "Default values:" << endl
         << "  > Document File: " << DEFAULT_DOC_FILE << endl
         << "  > Num. Clusters: " << DEFAULT_K << endl
         << "  > Num. Threads:  " << DEFAULT_THREADS << endl
         << "  > Mode:          " << "single thread (normal)" << endl
         << "    No vocabulary file (indices will be used instead)," << endl
         << "    using TXN scheme," << endl
         << "    displaying clustering results," << endl
         << "    optimization enabled." << endl;
    cout << "*To use max number of threads available, do not set t." << endl;
    cout << endl
         << "Example usage:" << endl
         << "  $ ./spkmeans -d ../TestData/news20 -k 50 --openmp" << endl;
    cout << endl;
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

        // find all documents in this partition and sum them together.
        float *sum = vec_zeros(data->wc);
        for(int j=0; j<(data->dc); j++) {
            if(data->p_asgns[j] == i) {
                for(int a=0; a<(data->docs[j].count); a++) {
                    float value = data->docs[j].words[a].value;
                    sum[data->docs[j].words[a].index] += value;
                }
            }
        }

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

        delete[] sum;
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
    unsigned int *k, unsigned int *num_threads, unsigned int *run_type,
    bool *use_scheme, bool *show_results, bool *auto_k, bool *optimize)
{
    // set defaults before proceeding to check arguments
    *doc_fname = DEFAULT_DOC_FILE;
    *vocab_fname = "";
    *k = DEFAULT_K;
    *num_threads = DEFAULT_THREADS;
    *run_type = RUN_NORMAL;
    *use_scheme = true;
    *show_results = true;
    *auto_k = false;
    *optimize = true;

    // check arguments: expected command as follows:
    // $ ./spkmeans -d docfile -w wordfile -k 2 -t 2 --galois
    for(int i=1; i<argc; i++) {
        string arg(argv[i]);

        // if flag is --help, return 1 to print usage instructions
        if(arg == "--help" || arg == "-help" || arg == "-h")
            return RETURN_HELP;
        // if flag is --version, return 2 to print version number
        else if(arg == "--version" || arg == "-version" || arg == "-V")
            return RETURN_VERSION;

        // if the flag was to run as galois or openmp, set the run type
        else if(arg == "--galois" || arg == "-galois")
            *run_type = RUN_GALOIS;
        else if(arg == "--openmp" || arg == "-openmp")
            *run_type = RUN_OPENMP;

        // also check if flag was set to disable weight normalization, squelch
        // displaying results, or disable optimizations
        else if(arg == "--noscheme" || arg == "-noscheme")
            *use_scheme = false;
        else if(arg == "--noresults" || arg == "-noresults")
            *show_results = false;
        else if(arg == "--noop" || arg == "-noop")
            *optimize = false;

        // or if K should be selected automatically
        else if(arg == "--autok" || arg == "-autok" ||
                arg == "--auto" || arg == "-auto")
            *auto_k = true;

        // otherwise, check the given flag value
        else {
            i++;
            if(i >= argc) {
                cout << "Warning: expected value after \"" << arg
                     << "\" argument. Continuing anyway." << endl;
                continue;
            }
            if(arg == "-d") // document file
                *doc_fname = string(argv[i]);
            else if(arg == "-w" || arg == "-v") // words file
                *vocab_fname = string(argv[i]);
            else if(arg == "-k") // size of k
                *k = atoi(argv[i]);
            else if(arg == "-t") // number of threads
                *num_threads = atoi(argv[i]);
            else { // otherwise, invalid input so print and decrement i again
                cout << "Unknown argument: \"" << arg
                     << "\". Use argument --help for more info." << endl;
                i--;
            }
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



// main: set up and start the clustering process.
int main(int argc, char **argv)
{
    // get file names, and set up k and number of threads
    string doc_fname, vocab_fname;
    unsigned int k, num_threads, run_type;
    bool use_scheme, show_results, auto_k, optimize;
    int retval = processArgs(argc, argv,
        &doc_fname, &vocab_fname, &k, &num_threads, &run_type,
        &use_scheme, &show_results, &auto_k, &optimize);
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

    // if K should be set automatically, compute that here if possible
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
#ifndef NO_GALOIS
    if(run_type == RUN_GALOIS) {
        // tell Galois the max thread count
        SPKMeansGalois spkm_galois(D, k, dc, wc, num_threads);
        if(!optimize)
            spkm_galois.disableOptimization();
        if(!use_scheme)
            spkm_galois.setScheme(SPKMeans::NO_SCHEME);
        cout << " [Galois: " << spkm_galois.getNumThreads()
             << " threads]." << endl;
        data = spkm_galois.runSPKMeans();
    }
#else
    if(run_type == RUN_GALOIS) {
        cout << endl << endl << "Error: GALOIS is not available."
             << "Please re-compile with Galois to use the \"--galois\" option."
             << endl << endl;
        return 0;
    }
#endif
    else if(run_type == RUN_OPENMP) {
        // tell OpenMP the max thread count
        SPKMeansOpenMP spkm_openmp(D, k, dc, wc, num_threads);
        if(!optimize)
            spkm_openmp.disableOptimization();
        if(!use_scheme)
            spkm_openmp.setScheme(SPKMeans::NO_SCHEME);
        cout << " [OpenMP: " << spkm_openmp.getNumThreads()
             << " threads]." << endl;
        data = spkm_openmp.runSPKMeans();
    }
    else {
        SPKMeans spkm(D, k, dc, wc);
        if(!optimize)
            spkm.disableOptimization();
        if(!use_scheme)
            spkm.setScheme(SPKMeans::NO_SCHEME);
        cout << " [single thread]." << endl;
        data = spkm.runSPKMeans();
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
