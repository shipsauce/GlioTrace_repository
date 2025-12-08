function A = count_a_raw(X)
% labelSeqs: cell array, each cell is a vector of labels (numeric or categorical)
%
% A(i,j) = number of times label i is followed by label j

    % Collect all unique labels in order
    allLabels = unique(vertcat(X{:}));

    % Init transition count matrix
    A = zeros(6, 6);

    % Count transitions for each sequence
    for s = 1:numel(X)
        seq = X{s};
        if numel(seq) < 2
            continue
        end

        % Count transitions seq(t) -> seq(t+1)
        for t = 1:numel(seq) - 1
            i = seq(t);
            j = seq(t+1);
            A(i, j) = A(i, j) + 1;
        end
    end
end
