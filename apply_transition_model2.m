function sequences = apply_transition_model2(sequences, A, B, pi)
% Run Viterbi algorithm on each observation sequence, handling NaN values,
% and updating the sequences with the inferred state sequences.
%
% Input parameters:
%   sequences - Matrix of observation sequences with NaN values (TxC)
%   A         - State transition matrix (NxN)
%   B         - Observation probability matrix (NxM)
%   pi        - Initial state distribution (1xN)
%
% Output parameters:
%   sequences - Updated matrix of sequences with inferred states

    [T, C] = size(sequences); % T: time steps, C: number of cells

    for c = 1:C
        % Extract the observations for the c-th cell, ignoring NaNs
        observations = sequences(:, c);
        validIdx = ~isnan(observations); % Find valid (non-NaN) entries
        validObservations = observations(validIdx);

        if ~isempty(validObservations)
            % Apply Viterbi algorithm only on valid observations
            inferredStates = viterbi(validObservations, A, B, pi);
            
            % Update the original sequences matrix with inferred states
            sequences(validIdx, c) = inferredStates;
        end
    end
end

