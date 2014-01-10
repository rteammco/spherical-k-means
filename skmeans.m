%%%%%%%%%% START HERE %%%%%%%%%%
disp('------------------------------------------------------------------');
disp('>>> Starting... <<<')



%%%%%%%%%% READ FILES %%%%%%%%%%
%files = {'./doc1.txt', './doc2.txt', './doc3.txt'};
files = {'./ds1.txt', './ds2.txt', './ds3.txt', 'ds4.txt', 'doc1.txt', 'doc2.txt', 'doc3.txt'};

wc = 0;
dc = 0;

vectors = {};
words = {};

for i = 1:length(files)
    doc = files{i};
    fh = fopen(files{i});
    if ~fh
        disp('>>> Could not open file. <<<');
        continue;
    end
    
    dc = dc + 1;
    doc_vector = [];
    
    while ~feof(fh)
        raw = fscanf(fh, '%s', 1);
        word = lower(regexprep(raw, '[^\w]', ''));
        
        if ~isempty(word)
            
            % TODO - eliminate common words
            
            % check if word already exists in the list
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
            
            % extend document vector as necessary (pad with 0s)
            %len = length(doc_vector);
            %while len < index
            %    len = len + 1;
            %    doc_vector{len} = 0;
            %end
            if length(doc_vector) < index
                doc_vector(index) = 0;
            end
            
            doc_vector(index) = doc_vector(index) + 1;
            
        end
        
    end
    
    fclose(fh);
    
    vectors{dc} = doc_vector;
    
end



%for i=1:length(vectors)
%   disp(vectors{i}); 
%end
%disp('---------------------------');

% set up values for all document vectors (txn scheme)
%
%txn scheme:
%				(wc) d = number of words in vector
%				f_ji = count of word j in document i
%				t_ji = f_ji
%				g_j = 1
%                   i fixed
%				s_i = sqrt( sum j=(1, d) { pow((t_ji * g_j), 2) } )
%
% x_i[j] = t_ji * g_j * s_i
for i = 1:length(vectors)
    
    
    % pad the vector to fit max length
    %len = length(vectors{i});
    %while len < wc
    %    len = len + 1;
    %    vectors{i}{len} = 0;
    %end
    if length(vectors{i}) < wc
        vectors{i}(wc) = 0;
    end
    
    
    % calculate s_i for this vector
    sum_i = 0;
    for k = 1:wc
        t_ki = vectors{i}(k);
        %disp(t_ki);
        q_k = 1;
        sum_i = sum_i + power((t_ki * q_k), 2);
    end
    %disp(sum_i);
    s_i = sqrt(sum_i);
    %disp(['>>> ' num2str(s_i)]);
    
    
    % calculate vector values
    for j = 1:wc
        t_ji = vectors{i}(j); % count of word j in document i
        q_j = 1; % constant for txn
        vectors{i}(j) = t_ji * q_j * s_i;
    end
    
    
end



%for i=1:length(vectors)
%   disp(vectors{i}); 
%end

disp(['Word Count: ' num2str(wc)]);
disp(['Document Count: ' num2str(dc)]);



%%%%%%% FIRST PARTITION %%%%%%%%

% partition into k disjoint clusters
k = 3;

% start with an arbitrary partitioning
split = floor(dc / k);
disp(['Split = ' num2str(split)]);

partitions = cell(1, k); % list of k partitions

base = 1;
for i = 1:k % for each partition
    top = base + split - 1;
    
    if i == k
        top = dc;
    end
    
    p_size = top - base + 1;
    partitions{i} = cell(1, p_size);
    
    disp(['Base = ' num2str(base)]);
    for j = 1:p_size % for each element in this partition
        v_loc = j + base - 1;
        %disp(v_loc);
        partitions{i}{j} = vectors{v_loc};
    end
    
    base = base + split;
end



% for all partitions, compute the concept vector of said partition
concepts = cell(1, k);
for i = 1:k
    n_i = length(partitions{i}); % number of document vectors in partition
    
    % get sum of all vectors in partition i
    sum_vector = zeros(1, wc);
    for j = 1:n_i
        for a = 1:wc
             sum_vector(a) = sum_vector(a) + partitions{i}{j}(a);
        end
    end
    %disp(sum_vector);
    
    % compute mean vector for partition and its norm
    m_i = (1 / n_i) * sum_vector;
    %disp(m_i);
    m_i_norm = norm(m_i);
    %disp(m_i_norm);
    
    % compute concept vector for partition using mean vector:
    c_i = m_i / m_i_norm;
    concepts{i} = c_i;
    %disp(c_i);
    
    %disp('---------');
end



threshold = 100;
Q = 0;
last_Q = -2 * threshold;
disp('-------STARTING ALGORITHM-------');
disp(['1st Q = ' num2str(Q)]);
disp(['1st last_Q = ' num2str(last_Q)]);
for i = 1:k
    n_i = length(partitions{i}); % num. doc vectors in partition
    disp(['1st size of partition ' num2str(i) ' = ' num2str(n_i)]);
end
disp('--------------------------------');

loopCounter = 0;
while abs(Q - last_Q) > threshold

    loopCounter = loopCounter + 1;
    disp(['LOOP #' num2str(loopCounter)]);
    
    % for each document vector, find closest concept vector of cosine
    %   similarity (dot product)
    new_partitions = cell(1, k);
    for i = 1:dc
        
        closest = 1;
        closest_val = -1;
        for j = 1:length(concepts)
            dot_p = dot(vectors{i}, concepts{j});
            norm_v = norm(vectors{i});
            norm_c = norm(concepts{j});
            
            similarity = dot_p / (norm_v * norm_c);
            if closest_val == -1 || similarity < closest_val
                closest_val = similarity;
                closest = j;
            end
        end
        
        % put vector into partition of the closest similarity
        new_partitions{closest}{end+1} = vectors{i};
        
    end
    
    
    % compute new concept vectors for the new_partitions set
    new_concepts = cell(1, k);
    for i = 1:k
        n_i = length(new_partitions{i}); % num. doc vectors in partition

        % get sum of all vectors in new partition i
        sum_vector = zeros(1, wc);
        for j = 1:n_i
            for a = 1:wc
                 sum_vector(a) = sum_vector(a) + new_partitions{i}{j}(a);
            end
        end

        % compute mean vector for partition and its norm
        m_i = (1 / n_i) * sum_vector;
        m_i_norm = norm(m_i);

        % compute concept vector for partition using mean vector:
        c_i = m_i / m_i_norm;
        new_concepts{i} = c_i;
        %disp(c_i);

        disp(['   Size of partition ' num2str(i) ' = ' num2str(n_i)]);
    end
    
    
    % compute new quality of partitioning
    last_Q = Q;
    Q = 0;
    for i = 1:k
        len = length(new_partitions{i});
        for j = 1:len % for each doc vector in partition i
            % sum += dot product of doc vector and partition i concept
            Q = Q + dot(new_partitions{i}{j}, new_concepts{i});
        end
    end
    disp(['   last_Q = ' num2str(last_Q)]);
    disp(['   Q = ' num2str(Q)]);
    
    % current lists become the new lists
    concepts = new_concepts;
    
    disp(['   > diff = ' num2str(abs(Q - last_Q))]);
    
end


disp('Partitions computed. End of code.');

partitions = new_partitions;


%%% LOOP %%%

% for each document vector, find concept vector of closest cosine
% similarity (??? dot product) TODO - what does closest mean?
%
% put vector into partition of the closest concept vector (random if tied);
% these will be new partitions at time t+1.

% compute new concept vectors

% if (stopping criterion met):
%   final_partitions = partitions[t+1];
%   break; % loop over, partition found!
%
% %"stopping criterion may be difference in Q at t and t+1 <= threshold."
% %TODO - how to calculate Q?

% t++;

%%% END LOOP %%%


% Document Vector = list of words and counts for a particular document
% 
% Mean Vector (in a cluster): m = (1 / #_doc_vecs_in_i) * (#_all_vecs_in_i)
% Concept Vector: c = m / ||m||
%
% 2) for all document vectors x, find closest concept vector in cosine
% similarity.
%
% 3) 
%

%
%(1) Start with arbitrary partitioning 0: {ℼ(0)j}kj=1.
%		{c(0)j}kj=1 denotes concept vectors associated with given partitioning 0.
%		Set t = 0 (index of iteration).
%	(2) ∀xi (1 to n), find concept vector closest in cosine similarity to xi.
%		Compute new partitioning {ℼ(t+1)j}kj=1 induced by old concept vectors {c(t)j}kj=1.
%		“ℼ(t+1)j is the set of all document vectors closest to the concept vector c(t)j.”
%		- if doc vector closest to more than one, assign it randomly to one of the clusters.
%	(3) Compute new concept vectors for (t+1):
%		c(t+1)j = m(t+1)j / || m(t+1)j ||, 1 ≤ j ≤ k.						( 8 )
%		m(t+1)j = centroid or mean of document vectors in cluster ℼ(t+1)j.
%	(4) 	If stopping criterion met: set ℼ*j = ℼ(t+1)j and c*j = c(t+1)j (for all 1 ≤ j ≤ k). EXIT.
%		Else: t++; GOTO step_2.
%		stopping criterion may be difference in Q at t and t+1 <= threshold.
%

















































































