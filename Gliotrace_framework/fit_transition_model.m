function [A, pi] = fit_transition_model(sequences, A, B, pi, maxIter, tol, mu, lambda)
% Baum-Welch algorithm for multiple observation sequences with known B
%
% Input parameters:
%   sequences - Cell array of observation sequences
%   A         - Initial state transition matrix (NxN)
%   B         - Known observation probability matrix (NxM)
%   pi        - Initial state distribution (1xN)
%   maxIter   - Maximum number of iterations
%   tol       - Tolerance for convergence
%   mu        - Strength of the diagonal prior
%   lambda    - L1-regularization parameter for sparsity
%
% Output parameters:
%   A         - Updated state transition matrix
%   pi        - Updated initial state distribution

    N = size(A, 1); % Number of states
    M = size(B, 2); % Number of possible observations

    prevLogLikelihood = -inf;

    for iter = 1:maxIter
        A_num = zeros(N, N);
        A_denom = zeros(N, 1);
        pi_num = zeros(1, N);
        logLikelihood = 0;

        for seqIdx = 1:length(sequences)
            observations = sequences{seqIdx};
            T = length(observations);

            % E-step: Compute alpha, beta, gamma, and xi
            alpha = zeros(T, N);
            beta = zeros(T, N);
            gamma = zeros(T, N);
            xi = zeros(T-1, N, N);

            % Forward algorithm
            alpha(1, :) = pi .* B(:, observations(1))';
            alpha(1, :) = alpha(1, :) / sum(alpha(1, :));
            for t = 2:T
                alpha(t, :) = (alpha(t-1, :) * A) .* B(:, observations(t))';
                alpha(t, :) = alpha(t, :) / sum(alpha(t, :));
            end

            % Backward algorithm
            beta(T, :) = 1;
            for t = T-1:-1:1
                beta(t, :) = (A * (B(:, observations(t+1)) .* beta(t+1, :)'))';
                beta(t, :) = beta(t, :) / sum(beta(t, :));
            end

            % Compute gamma
            for t = 1:T
                gamma(t, :) = alpha(t, :) .* beta(t, :);
                gamma(t, :) = gamma(t, :) / sum(gamma(t, :));
            end

            % Compute xi
            for t = 1:T-1
                denom = (alpha(t, :) * A) .* B(:, observations(t+1))' * beta(t+1, :)';
                for i = 1:N
                    for j = 1:N
                        xi(t, i, j) = alpha(t, i) * A(i, j) * B(j, observations(t+1)) * beta(t+1, j) / denom;
                    end
                end
            end

            % Accumulate statistics
            pi_num = pi_num + gamma(1, :);
            for i = 1:N
                A_denom(i) = A_denom(i) + sum(gamma(1:T-1, i));
                for j = 1:N
                    A_num(i, j) = A_num(i, j) + sum(xi(:, i, j));
                end
            end

            % Compute log-likelihood for convergence check
            logLikelihood = logLikelihood + sum(log(sum(alpha, 2)));
        end

        % M-step: Update pi and A with regularization
        pi = pi_num / sum(pi_num);

        for i = 1:N
            for j = 1:N
                num = A_num(i, j) + mu * (i == j) - lambda * (i ~= j);
                denom = A_denom(i) + mu + lambda;
                A(i, j) = max(0, num / denom); % Ensure non-negativity
            end
            A(i, :) = A(i, :) / sum(A(i, :)); % Re-normalize to ensure valid probabilities
        end

        % Check for convergence
        if abs(logLikelihood - prevLogLikelihood) < tol
            break;
        end
        prevLogLikelihood = logLikelihood;
    end
end