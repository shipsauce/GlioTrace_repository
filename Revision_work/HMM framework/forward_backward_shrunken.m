function [gamma, xi, log_lik] = forward_backward_shrunken( ...
        X, trajectory, A_raw, mu_shrunk, sigma2_iso, active_mask, glm_models, pi_vec, eps)
% @authors: André Lasses Armatowski, Linnea Hallin 
% @date: 09112025

if nargin < 9
    eps = 1e-12;
end

T = size(trajectory, 1);
K = size(mu_shrunk, 1);

log_alpha = zeros(T, K);
log_beta  = zeros(T, K);
log_xi    = zeros(T-1, K, K);
log_gamma = zeros(T, K);

% Log emissions: (T, K)
log_em = log_emission_shrunken_iso(X, mu_shrunk, sigma2_iso, active_mask, eps);

% Log transitions: (T, K, K)
log_A = evaluate_models_on_trajectories(glm_models, A_raw, trajectory, K, eps);

% Log initial distribution
log_pi = log(max(pi_vec(:), eps));
log_pi = log_pi - logsumexp_vec(log_pi);

% Forward recursion
log_alpha(1, :) = log_pi.' + log_em(1, :);

for t = 2:T
    % log_alpha(t, j) = log_em(t, j) + logsumexp_i( log_alpha(t-1, i) + logA(t-1,i,j) )
    prev = reshape(log_alpha(t-1, :), [K 1]);
    trans = squeeze(log_A(t-1, :, :));                 % (K, K)
    tmp = prev + trans;                                % (K, K)
    log_alpha(t, :) = log_em(t, :) + logsumexp_mat(tmp, 1);
end

% Log likelihood
log_lik = logsumexp_vec(log_alpha(T, :).');

% Backward recursion
log_beta(T, :) = 0;

for t = T-1:-1:1
    trans  = squeeze(log_A(t, :, :));                  % (K, K)
    emit   = log_em(t+1, :);                           % (1, K)
    nxt    = log_beta(t+1, :);                         % (1, K)
    tmp    = trans + repmat(emit + nxt, K, 1);         % (K, K)
    log_beta(t, :) = logsumexp_mat(tmp, 2).';
end

% Gamma
log_gamma = log_alpha + log_beta;
ls = logsumexp_mat(log_gamma, 2);
log_gamma = log_gamma - ls;
gamma = exp(log_gamma);

% Xi
for t = 1:T-1
    la = reshape(log_alpha(t, :), [K 1]);
    trans = squeeze(log_A(t, :, :));                  % (K, K)
    emit = log_em(t+1, :);
    lb = log_beta(t+1, :);
    tmp = la + trans + repmat(emit + lb, K, 1);       % (K, K)
    tmp = tmp - logsumexp_vec(tmp(:));
    log_xi(t, :, :) = tmp;
end

xi = exp(log_xi);
end
