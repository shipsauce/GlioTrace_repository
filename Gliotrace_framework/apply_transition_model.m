function viterbiPaths = apply_transition_model(sequences, A, B, pi)
% Run Viterbi algorithm on each observation sequence
%
% Input parameters:
%   sequences - Cell array of observation sequences
%   A         - State transition matrix (NxN)
%   B         - Observation probability matrix (NxM)
%   pi        - Initial state distribution (1xN)
%
% Output parameters:
%   viterbiPaths - Cell array of inferred state sequences for each observation sequence

    numSeq = length(sequences);
    viterbiPaths = cell(numSeq, 1);

    for i = 1:numSeq
        observations = sequences{i};
        viterbiPaths{i} = viterbi(observations, A, B, pi);
    end
end

