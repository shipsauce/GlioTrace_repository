function [pi, mu_shrunk, sigma2_iso, active_mask, glm_models, A_global] = ...
    HMM_glm_shrunken(trajectories, embeddings, hard_labels, pi, glm_models, mu_shrunk, sigma2_iso, active_mask, ...
                     K, max_iter, keep_pct, glm_iters, penalty, eps_conv)
% @authors: André Lasses Armatowski, Linnea Hallin 
% @date: 09112025

    if nargin < 9, K = 6; end
    if nargin < 10, max_iter = 10; end
    if nargin < 11, keep_pct = 0.1; end
    if nargin < 12, glm_iters = 500; end
    if nargin < 13, penalty = 1; end
    if nargin < 14, eps_conv = 1e-6; end

    prev_lik = -inf;
    N_seq = numel(trajectories);
    
    A_raw = count_a_raw(hard_labels);

    for it = 1:max_iter
        gammas = cell(1, N_seq);
        xis = cell(1, N_seq);
        full_lik = 0;

        % ---- E-step ----
        for k = 1:N_seq
            X = embeddings{k};
            F = trajectories{k};

            [gamma, xi, log_lik] = forward_backward_shrunken(X, F, A_raw, mu_shrunk, sigma2_iso, active_mask, glm_models, pi);

            gammas{k} = gamma;
            xis{k} = xi;
            full_lik = full_lik + log_lik;
        end

        % Early convergence check
        if abs(full_lik - prev_lik) < eps_conv
            disp('Early convergence');
            return;
        end
        prev_lik = full_lik;

        % ---- M-step ----
        % Update pi: average initial posterior across sequences
        pi_stack = cell2mat(cellfun(@(g) g(1,:), gammas', 'UniformOutput', false));
        pi = mean(pi_stack, 1);
        pi = max(pi, 1e-12);
        pi = pi / sum(pi);

        A_global = calc_A_global(xis);

        % Update emissions: shrunken centroids, shared variance, active mask
        [mu_shrunk, sigma2_iso, active_mask] = update_emission_shrunken(embeddings, gammas, keep_pct);

        % Update transitions: multinomial logistic per "from" state
        glm_models = update_transitions(trajectories, xis, glm_iters, ...
                                                    penalty, glm_models);
        it
    end
end
