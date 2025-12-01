function A = count_a_raw(X)
% labelSeqs: cell array, each cell is a vector of labels (numeric or categorical)
%
% A(i,j) = number of times label i is followed by label j

    % Collect all unique labels in order
    allLabels = unique(vertcat(X{:}));
    L = numel(allLabels);

    % Map labels to 1..L
    % Works for numeric or categorical
    [~, ~, idxGlobal] = unique(allLabels);

    % Build lookup from label value to index
    label2idx = containers.Map(num2cell(allLabels), num2cell(idxGlobal));

    % Init transition count matrix
    A = zeros(L, L);

    % Count transitions for each sequence
    for s = 1:numel(X)
        seq = X{s};
        if numel(seq) < 2
            continue
        end

        % Map sequence to indices
        idxSeq = arrayfun(@(x) label2idx(x), seq);

        % Count transitions idxSeq(t) -> idxSeq(t+1)
        for t = 1:numel(idxSeq) - 1
            i = idxSeq(t);
            j = idxSeq(t+1);
            A(i, j) = A(i, j) + 1;
        end
    end
end
