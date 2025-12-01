function [mu_shrunk, sigma2_iso, active_mask] = update_emission_shrunken(Xs, gammas, keep_pct, eps)
    if nargin < 3
        keep_pct = 0.1;
    end
    if nargin < 4
        eps = 1e-12;
    end

    K = size(gammas{1}, 2);
    D = size(Xs{1}, 2);

    % Stack for vectorized operations
    X_all = vertcat(Xs{:});       % (T_total, D)
    G_all = vertcat(gammas{:});   % (T_total, K)

    % Raw class means
    Nks = sum(G_all, 1)';         % (K, 1)
    means_raw = (G_all' * X_all); %weighted sum of the data points for each class
    means_raw = means_raw ./ (Nks + eps);  % (K, D)
    Ntot = sum(Nks) + eps;

    % Global mean
    mu_global = (Nks .* means_raw)' * ones(K,1) / Ntot;
    mu_global = mu_global(:)';  % (1, D)

    % Compute squared differences
    % diff: (T_total, K, D)
    T_total = size(X_all,1);
    diff = reshape(X_all, [T_total,1,D]) - reshape(means_raw, [1,K,D]);

    % Pooled variance per feature
    G_all_exp = reshape(G_all, [T_total,K,1]);
    s2 = sum(sum(G_all_exp .* diff.^2, 1), 2) / Ntot;
    s2 = reshape(s2, [1,D]);

    % Regularize and sqrt
    s = sqrt(max(s2, eps));

    % Standardized class mean deviations
    d = (means_raw - mu_global) ./ s;  % (K, D)

    % Feature importance
    feat_importance = mean(abs(d), 1);

    % Threshold top features
    thresh = prctile(feat_importance, 100*(1-keep_pct));
    active_mask = feat_importance >= thresh;

    % Shrink
    d_tilde = d .* active_mask;

    % Reconstruct shrunk means
    mu_shrunk = mu_global + d_tilde .* s;

    % Active feature count
    D_eff = sum(active_mask);

    % Shared isotropic variance over active dims
    diff_active = reshape(X_all(:,active_mask), [T_total,1,D_eff]) - reshape(mu_shrunk(:,active_mask), [1,K,D_eff]);
    num = sum(sum(sum(G_all_exp .* diff_active.^2,1),2),3);
    sigma2_iso = num / (D_eff * Ntot);
end
