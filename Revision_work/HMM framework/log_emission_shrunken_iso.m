function log_prob = log_emission_shrunken_iso(X, mu_shrunk, sigma2, active_mask, eps)
% Vectorized log-emission for shrunken-centroid isotropic Gaussian.
% X: (T, D)
% mu_shrunk: (K, D)
% sigma2: scalar
% active_mask: logical vector (D,) or []
% Returns: (T, K)
%
%
% @authors: André Lasses Armatowski, Linnea Hallin 
% @date: 09112025

if nargin < 5
    eps = 1e-6;
end

if ~isempty(active_mask)
    X_ = X(:, active_mask);
    M_ = mu_shrunk(:, active_mask);
else
    X_ = X;
    M_ = mu_shrunk;
end

[T, D_eff] = size(X_);
K = size(M_, 1);

% Expand for broadcasting: X_ -> (T, 1, D_eff), M_ -> (1, K, D_eff)
X_exp = reshape(X_, [T 1 D_eff]);
M_exp = reshape(M_, [1 K D_eff]);

diff = X_exp - M_exp;                     % (T, K, D_eff)
dist2 = sum(diff.^2, 3);                  % (T, K)

log_norm = -0.5 * D_eff * log(2 * pi * (sigma2 + eps));
log_prob = log_norm - 0.5 * dist2 / (sigma2 + eps);

end
