% WARNING:
% This is an inefficient implementation. Run once, and continue
% using the data set.


function [ D, words ] = readDocsFromDir( path, fillerFile )
% readDocs: Reads documents from the given directory and returns them as
%   a sparse matrix (D) whose columns represent documents.
% PARAMETERS:
%   path - directory which contains all documents to read (text files).
%   fillerFile - path/filename of a file containing filler words.
% RETURNS:
%   D - sparse matrix containing all document vectors.
%   words - a list of all words (indexed identically to D).


    words = {};     % list of words
    D = zeros(1,0); % document matrix
    first = 1;      % first loop tracker
    wc = 0;         % word count tracker
    
    
    % Read file of "filler" words and aggregate filler words list:
    filler_words = {};
    if exist(fillerFile)
        filler_fh = fopen(fillerFile);
        if filler_fh ~= -1
            while ~feof(filler_fh)
                filler_words{end+1} = fscanf(filler_fh, '%s', 1);
            end
            fclose(filler_fh);
        end
    end
    
    
    % Read each file and attempt to parse it into the document matrix:
    files = dir(path);
    num_files = length(files);
    
    % Fix out-of-memory issues
    if num_files > 2000
        num_files = 2000;
    end
    
    % Shuffle the file list so documents aren't already pseudo-clustered:
    files = files(randperm(num_files));
    
    for i = 1:num_files
        
        disp(['Reading file ' num2str(i) ' of ' num2str(num_files) ...
              ': ' files(i).name]);
        
        % Skip directories (only interested in files):
        if files(i).isdir
            continue
        end
        
        % Try to open the file (skip if can't open):
        fname = files(i).name;
        fh = fopen([path '/' fname]);
        if ~fh
            disp(['READ ERROR: Could not open file: ' fname]);
            continue;
        end
        
        
        % Create the document vector (pre-allocate zeros):
        [lenD, widD] = size(D);
        doc_vector = zeros(1, lenD);
        

        % Read the file word-by-word, and add words to the vector:
        while ~feof(fh)
            raw = fscanf(fh, '%s', 1);
            
            % Strip non-alphanumeric characters and convert to lowercase
            word = lower(regexprep(raw, '[^\w]', ''));

            % Process word if anything is left after cleaning
            if ~isempty(word)

                % Eliminate common "junk" or "filler" words
                [isFiller, ~] = ismember(word, filler_words);
                if isFiller
                    continue;
                end

                % Find index where this word goes; if the word was not
                %   seen before, expand the list and add it to the end
                index = 0;
                for i = 1:length(words)
                    if strcmp(words{i}, word)
                        index = i;
                        break
                    end
                end
                if index == 0
                    wc = wc + 1;
                    words{wc} = word;
                    index = wc;
                    doc_vector(index) = 0;
                end

                % Increment this word's count in this document vector
                doc_vector(index) = doc_vector(index) + 1;

            end

        end
        fclose(fh);
        
        
        if first
            D = doc_vector';
            first = 0;
        else
            % Pad vector or matrix with zeros to match sizes if necessary:
            lenV = length(doc_vector);
            if lenV < lenD
                doc_vector(lenD) = 0;
            elseif lenD < lenV
                D(lenV, widD) = 0;
            end

            % Add this vector the the document matrix:
            D = [D doc_vector'];
        end
        
    end


end

