% WARNING:
% This is an inefficient implementation. Run once, and continue
% using the data set.


function [ D, words ] = readDocsFromFile( docFile, vocabFile )
% readDocs: Reads documents from the given document file and vocabulary
%   file and returns them as a sparse matrix (D) whose columns represent
%   documents.
% PARAMETERS:
%   docFile - path to file containing the document list.
%   vocabFile - path to a vocabulary (dictionary) file containing a
%       list of all the words used in the documents.
% RETURNS:
%   D - sparse matrix containing all document vectors.
%   words - a list of all words (indexed identically to D).


    % Try to open the document file:
    if ~exist(docFile, 'file')
        disp(['ERROR: file not found: ' docFile]);
        return;
    end
    docFH = fopen(docFile);
    if ~docFH
        disp(['ERROR: could not open: ' docFile]);
        return;
    end
    
    % Try to open the vocabulary file:
    if ~exist(vocabFile, 'file')
        disp(['ERROR: file not found: ' vocabFile]);
        return;
    end
    vocabFH = fopen(vocabFile);
    if ~vocabFH
        disp(['ERROR: could not open: ' vocabFile]);
        return;
    end
    
    
    % Read metadata (num. words, num. docs from docFile)
    numDocs = fscanf(docFH, '%d', 1);
    wc = fscanf(docFH, '%d', 1);
    
    disp(['Docs: ' num2str(numDocs) '  WC: ' num2str(wc)]);

    
    % Read in the vocabulary words:
    words = cell(1, wc);
    curWord = 1;
    while ~feof(vocabFH)
        word = fscanf(vocabFH, '%s', 1);
        if ~isempty(word)
            words{curWord} = word;
            curWord = curWord + 1;
        end
    end
    fclose(vocabFH);
    
    
    % Read in the document matrix:
    D = zeros(wc, numDocs);
    fscanf(docFH, '%s', 1); % get rid of extra line
    while ~feof(docFH)
        docID = fscanf(docFH, '%d', 1);
        wordID = fscanf(docFH, '%d', 1);
        count = fscanf(docFH, '%d', 1);
        D(wordID, docID) = count;
    end
    fclose(docFH);


end

