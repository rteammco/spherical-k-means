/* File: reader.cpp
 *
 * Definitions of the document and vocabulary file reading functions.
 */

#include "reader.h"

#include <fstream>
#include <sstream>
#include <string>
#include <string.h>

using namespace std;



// Read the document data file into a spare matrix (2D array) format).
// This function assumes that the given file name is valid.
// Returns the 2D array (sparse matrix) - columns are document vectors.
float** readDocFile(const char *fname, int &dc, int &wc)
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

    infile.close();
    return mat;
}



// Read the word data into a list. Words are just organized one word per line.
// Returns a list of strings (char pointers), or a null pointer if the given
// file name does not exist.
char** readWordsFile(const char *fname, int wc)
{
    ifstream infile(fname);

    // check that the file exists - if not, return a null pointer
    if(!infile.good()) {
        infile.close();
        return 0;
    }

    // read in the list of words
    char **words = new char*[wc];
    int i = 0;
    string line;
    while(getline(infile, line)) {
        // if we're out of word space and the file has more, stop
        if(i >= wc)
            break;

        // copy the word into memory
        char *word = new char[line.size() + 1];
        word[line.size()] = 0;
        memcpy(word, line.c_str(), line.size());
        words[i] = word;
        i++;
    }

    infile.close();
    return words;
}
