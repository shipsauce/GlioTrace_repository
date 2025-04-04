function path = viterbi(observations, A, B, pi)
% Viterbi algorithm to find the most likely sequence of hidden states
%
% Inputs:
%   observations - Sequence of observations (1xT)
%   A            - State transition matrix (NxN)
%   B            - Observation probability matrix (NxM)
%   pi           - Initial state distribution (1xN)
%
% Outputs:
%   path - Most likely sequence of hidden states (1xT)

    N = size(A, 1);
    T = length(observations);

    delta = zeros(T, N);
    psi = zeros(T, N);

    % Initialization
    delta(1, :) = pi .* B(:, observations(1))';
    psi(1, :) = 0;

    % Recursion
    for t = 2:T
        for j = 1:N
            [delta(t, j), psi(t, j)] = max(delta(t-1, :) .* A(:, j)' * B(j, observations(t)));
        end
    end

    % Termination
    [~, path(T)] = max(delta(T, :));

    % Path backtracking
    for t = T-1:-1:1
        path(t) = psi(t+1, path(t+1));
    end
end