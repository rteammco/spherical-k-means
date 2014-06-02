/* File: spkmeans.cpp
 *
 * A parallel implementation of the Spherical K-Means algorithm using the
 * Galois library (URL_HERE).
 */


#include <fstream>
#include <iostream>
#include <math.h>
#include <string>

#include "Galois/Galois.h"

#include "vectors.h"


#define DEFAULT_K 2;
#define DEFAULT_THREADS 2;


using namespace std;



// Runs the spherical k-means algorithm on the given sparse matrix D and
// clusters the data into k partitions.
void runSPKMeans(float **doc_matrix, unsigned int k, int dc, int wc)
{
    // TODO - make these parameters
    float qThresh = 0;

    // normalize each of the document vectors
    for(int i=0; i<dc; i++)
        vec_normalize(doc_matrix[i], wc);

    cout << "Running spherical K-means." << endl;
}



/* Read the document data file into a spare matrix (2D array) format).
 * This function assumes that the given file name is valid.
 * Returns the 2D array (sparse matrix) - columns are document vectors.
 * FILE FORMAT:
 *  <top of file>
 *      number of documents
 *      number of unique words
 *      number of non-zero words in the collection
 *      docID wordID count
 *      docID wordID count
 *      ...
 *  <end of file>
 */
float** loadFile(const char *fname, int &dc, int &wc)
{
    ifstream infile(fname);
    
    // get the number of documents and words in the data set
    int nzwc; // non-zero word count (ignored)
    infile >> dc >> wc >> nzwc;

    // set up the matrix and initialize it to all zeros
    float **mat = new float*[dc];
    for(int i=0; i<dc; i++) {
        mat[i] = new float[wc];
        memset(mat[i], 0, wc * sizeof(*mat[i]));
    }

    // populate the matrix from the data file
    string line;
    while(getline(infile, line)) {
        istringstream iss(line);
        int doc_id, word_id, count;
        if(!(iss >> doc_id >> word_id >> count))
            continue;
        mat[doc_id-1][word_id-1] = count;
    }

    return mat;
}



// Takes argc and argv from program input and parses the parameters to set
// values for k (number of clusters) and num_threads (the maximum number
// of threads for Galois to use).
void processArgs(int argc, char **argv,
    string *fname, unsigned int *k, unsigned int *num_threads)
{
    // set up file name
    if(argc >= 2)
        *fname = string(argv[1]);
    else
        *fname = "data";

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
}



// main: set up Galois and start the clustering process.
int main(int argc, char **argv)
{
    // get file name, and set up k and number of threads
    string fname;
    unsigned int k, num_threads;
    processArgs(argc, argv, &fname, &k, &num_threads);

    // tell Galois the max thread count
    Galois::setActiveThreads(num_threads);
    num_threads = Galois::getActiveThreads();
    cout << "Running SPK Means on " << fname << " with k=" << k
         << " (" << num_threads << " threads)." << endl;

    // set up the sparse matrix
    int dc, wc;
    float **D = loadFile(fname.c_str(), dc, wc);

    // run spherical k-means on the given sparse matrix D
    runSPKMeans(D, k, dc, wc);

    return 0;
}
