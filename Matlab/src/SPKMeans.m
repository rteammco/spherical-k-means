function [ P, cVs ] = SPKMeans( D, k )
% SPKMeans: Clusters the documents in given document matrix D using the
%   Spherical K-Means Algorithm.
% PARAMETERS:
%   D - the document vector matrix to be clustered.
%   k - number of clusters to partition documents into.
% RETURNS:
%   P - the set of partitions created by the algorithm.
%   cVs - the concept vectors of the returned partition.


    tic;
    
    
    % Set a threshold for change in quality stop condition.
    qThresh = 0.001;
    

    [~, num_docs] = size(D);
    for i = 1:num_docs 
        D(:,i) = D(:,i)/norm(D(:,i));
    end
    
    % Apply TNX or TFN scheme:
    D = txnScheme(D);
    
    
    
    % Create first partitioning (P):
    [wc, dc] = size(D);
    split = floor(dc / k);
    disp(['Split = ' num2str(split)]);

    P = cell(1, k);

    base = 1;
    for i = 1:k % for each partition
        top = base + split - 1;

        if i == k
            top = dc;
        end

        p_size = top - base + 1;
        %P{i} = cell(1, p_size);
        P{i} = zeros(wc, 0);

        %disp(['Base = ' num2str(base)]);
        for j = 1:p_size % for each element in this partition
            v_loc = j + base - 1;
            P{i} = [P{i} D(:, v_loc)];
            %P{i}{j} = D(:, v_loc);
        end

        base = base + split;
        %disp(P{i});
    end
    
    
    
    % Compute concepts for each partition:
    cVs = getConcepts(P, k);
    
    
    % Get initial partition quality:
    q = 0;
    for i = 1:k
        q = q + dot(sum(P{i}, 2), cVs{i});
    end
    dQ = qThresh * 10;
    t = 0;
    % Loop as long as the change in quality is significant:
    while dQ > qThresh && t < 200
        
        % Empty the partitions:
        for i=1:k
            P{i} = zeros(wc, 0);
        end
        
        
        % For each document vector, find closest concept vector of cosine
        %   similarity (smallest dot product):
        for i = 1:dc
            
            closest = 1;
            %cVal = dot(D(:,i), cVs{1}) / (norm(D(:,i)) * norm(cVs{1}));
            cVal = D(:,i)' * cVs{1};
            
            for j=2:k
                % Compute the similarity
                %s = dot(D(:,i), cVs{j}) / (norm(D(:,i)) * norm(cVs{j}));
                s = D(:,i)' * cVs{j};
                
                % If this concept is more similar, choose it instead.
                %   Similarity is closer when s is larger.
                if s > cVal
                    cVal = s;
                    closest = j;
                end
            end
            
            % Assign the vector to the closest partition
            P{closest} = [P{closest} D(:,i)];
        end
        
        
        % Compute the new partition's concepts:
        cVs = getConcepts(P, k);
        
        
        % Compute the quality of the current partitioning:
        lastQ = q;
        q = 0;
        for i = 1:k
            if ~isempty(P{i})
                q = q + dot(sum(P{i}, 2), cVs{i});
            end
        end
        dQ = abs(q - lastQ);
        disp(['q=' num2str(q) ' dQ=' num2str(dQ)]);
        
        
        t = t + 1;
        
    end
    
    
    toc;
    
    
end



function [ concepts ] = getConcepts( partitions, k )
% getConcepts: returns a set of concept vectors, one for each of the
%   given partitions. Concept vectors are indexed in the same order
%   as the partitions given.
% PARAMETERS:
%   partitions - the set of partitions of the document vectors.
%   k - number of partitions.
% RETURNS:
%   concepts - the set of concept vectors for each partition.

    % For all partitions, compute the concept vector of said partition:
    concepts = cell(1, k);
    for i = 1:k
        
        % If partition is empty, make the concept vector a zero vector:
        if isempty(partitions{i})
            [wc, ~] = size(partitions{i});
            concepts{i} = zeros(wc, 1);
            
        else
            % NOTE - mult. by 1/len(P{i}) not in some implementations...
            % Compute mean vector for partition:
            mV_i = (1 / length(partitions{i})) * sum(partitions{i}, 2);

            % Compute concept vector for partition using mean vector:
            concepts{i} = mV_i / norm(mV_i);
        end
        
    end

end



function [ D ] = txnScheme( D )
% txnScheme: Helper function for the SPKMeans algorithm. This function
%   document matrix D and modifies it with the tnx scheme.
% PARAMETERS:
%   D - un-adjusted document matrix D.
% RETURNS:
%   D - adjusted document matrix D.
%
% txn/tfn scheme algorithm (for each document vector x_i):
%      d = number of words in vector
%	f_ji = count of word j in document i
%	t_ji = f_ji
%	 g_j = 1 [txn]
%        = log (n / d_j) [tfn: d_j = # docs containing word j]
%	 s_i = sqrt( sum j=(1, d) { pow((t_ji * g_j), 2) } )
%
%   x_i[j] = t_ji * g_j * s_i


    [~, num_docs] = size(D);
    
    for i = 1:num_docs

        % Calculate s_i for this vector:
        D(:,i) = D(:,i) * sqrt(sum(D(:,i).^2));
        
    % NOTE: for tfn, construct a vector for d_j (# of documents
    %   each word j appears in. Then, multiply:
    % D(:,i) = times(D(:,i), d_vector) * s_i
        
        %vec_i = D(:,i);
        %tempv = vec_i.^2;
        %sum_i = sum(tempv);
        %s_i = sqrt(sum_i);
        %D(:,i) = D(:,i) * s_i
        
        %sum_i = 0;
        %for j = 1:wc
        %    f_ji = D(j, i);
        %    t_ji = f_ji;
        %    g_j = 1;
        %    %sum_i = sum_i + power((t_ji * g_j), 2);
        %    p = (t_ji * g_j);
        %    sum_i = sum_i + (p * p);
        %end
        %s_i = sqrt(sum_i);

        % calculate vector values
        %for j = 1:wc
        %    t_ji = vectors{i}(j); % count of word j in document i
        %    q_j = 1; % constant for txn
        %    vectors{i}(j) = t_ji * q_j * s_i;
        %end

    end
    
end