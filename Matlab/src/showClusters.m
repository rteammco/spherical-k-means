function [] = showClusters( P, D, words, numToShow )
% showClusters: ...

    if nargin < 4
        numToShow = 10;
    end

    k = length(P);
    for i = 1:k
        disp(['Partition: ' num2str(i)]);
        pSum = sum(P{i}, 2);
        [vals, indices] = sort(pSum, 'descend');
        for j = 1:numToShow
            if vals(j) == 0
                break;
            else
                disp(['   ' words{indices(j)}]);
            end
        end
    end
end