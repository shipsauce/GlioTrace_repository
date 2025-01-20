function subsequences = extractSubsequences(sequence, identity)
% This function takes a label sequence and returns the index of one or more
% subsequences where the label agrees with the dominating label of the
% sequence.
%
% Input parameters:
% sequence - label sequence
% identity - most prevalent label in the sequence
%
% Output parameters:
% subsequences - cell array of indicies of subsequences whose label agrees
% with the identity of the sequence
%
% @date: 16012025

% Initialize an empty cell array to store valid subsequences
subsequences = {};

% Initialize variables to keep track of current subsequence
currentSubsequence = [];

% Loop through the sequence to find subsequences that match the identity
for i = 1:length(sequence)
    if sequence(i) == identity
        % If the current element matches the identity, add it to the current subsequence
        % currentSubsequence = [currentSubsequence, sequence(i)];
        currentSubsequence = [currentSubsequence, i];
    else
        % If the current element does not match, check if the current subsequence is valid
        if length(currentSubsequence) >= 3
            subsequences{end+1} = currentSubsequence;
        end
        % Reset the current subsequence
        currentSubsequence = [];
    end
end

% Final check in case the sequence ends with a valid subsequence
if length(currentSubsequence) >= 3
    subsequences{end+1} = currentSubsequence;
end
end