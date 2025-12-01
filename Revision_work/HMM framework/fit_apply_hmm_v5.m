function fit_apply_hmm_v5(tbl)

% Add empty columns to table
tbl.A_est = cell(height(tbl),1);
tbl.pi_est = cell(height(tbl),1);
K = 6; % no of states

% Shrunken centroid initiation
init_params = readtable("shrunken_centroids_init.txt");
mu_init = table2array(init_params(1:6,:)); % initial centroid coords
s_j = table2array(init_params(7,:)); % dimension-wise std
sigma2_iso = sqrt(mean(s_j.^2));
active_mask = ones(1, size(s_j,2));

pi_init = ones(1, K) / K; % Uniform initial state distribution

% HMM hyperparams
max_iter = 50;
keep_pct=1;
glm_iters = 200;
penalty = 0.1;

cellines = unique(tbl.HGCC);

% Loop through the cellines
for i=1:length(cellines)
    fprintf(['Fit HMM parameters for celline: ' cellines{i} '...\n'])
    hgcc = cellines{i};
    subtable = tbl(tbl.HGCC == string(hgcc),:);

    perturbations = unique(subtable.perturbation);

    % Loop through the perturbations
    for j=1:length(perturbations)
        pert = perturbations{j};
        subtable_2 = subtable(subtable.perturbation == string(pert),:);
        dosez = unique(subtable_2.dose);

        % Loop through the doses
        for k=1:length(dosez)
            dose_curr = dosez(k);
            subtable_3 = subtable_2(subtable_2.dose == dose_curr,:);
            sequences = {};
            propagated_labels = {};
            tme_features = {};
            idx = logical((tbl.HGCC == string(hgcc)) .* (tbl.perturbation == string(pert)) .* (tbl.dose == dose_curr));

            [embeddings, hard_labels, propagated_labels, tme_features, delta_ts, comb_features] = convert_data_for_hmm(subtable_3);

    
            [pi, mu_shrunk, sigma2_iso, active_mask, glm_models, A_global] = HMM_glm_shrunken(tme_features, ...
                                                                                    embeddings, ...
                                                                                    hard_labels,...
                                                                                    pi_init, ...
                                                                                    [], ...
                                                                                    mu_init, ...
                                                                                    sigma2_iso, ...
                                                                                    active_mask, ...
                                                                                    K, ...
                                                                                    max_iter, ...
                                                                                    keep_pct, ...
                                                                                    glm_iters, ...
                                                                                    penalty);
    
            tbl.A_est(idx) = {A_t};
            tbl.pi_est(idx) = {pi_est};

        end
    end
end

end