% HMM with Shrunken-Centroid Emissions
% K latent states
% X{i}(t,:) = feature vector for cell i, time t
% E{i}(t,:) = embedding vector (observation)
% Transition probabilities depend on X via softmax
% Emission probability ∝ exp(–½ * squared standardized distance to shrunken centroid)

load("U3180_ctrl_tbl.mat")

% PARAMETERS & HYPERPARAMETERS
K = 6; 
tol = 1e-4; 
maxIter = 50; 
lambda = 1e-3; 
Delta = 0.5; 
feature_dependent_transitions = 0; % Toggle on/off

% LOAD / PREPARE DATA
[emissions, hard_labels, propagated_labels, tme_features] = convert_data_for_hmm(U3180_ctrl_tbl(U3180_ctrl_tbl.perturbation == "control",:)); % funktionen kan returnera för alla komb
E = emissions;

%% TEST ON SMALL SUBSET OF EMBEDDINGS
E = emissions(1:1000,:);
propagated_labels = propagated_labels(1:1000,:);
tme_features = tme_features(1:1000,:);
hard_labels = hard_labels(1:1000,:);
Ncells = numel(E);

if(feature_dependent_transitions)
    X = tme_features;
else
    for i=1:Ncells
        Ti = size(tme_features{i},1);
        X{i} = ones(Ti,1);
    end
    X = X';
end

d = size(E{1},2); % embedding size
p = size(X{1},2); % tme feature size

% INITIALIZATION
pi_k = ones(1,K) / K; % Equal state dist at start

% Initialize mu and MAD for shrunken centroid emissions
% 1) filter away dimensions with a high fraction of zeros
allE = cat(1, E{:}); 
zero_idx = allE==0; % find zeros in the matrix
frac_zero = mean(zero_idx,1); % compute the fraction of zeros for each dimension
keep = frac_zero < 0.98;    % threshold
X0 = allE(:, keep); % remove dimensions below threshold

% 2) robust center-scale on kept dims
med = median(X0,1);
madv = mad(X0,1); % compute median absolute deviation across dims
madv = max(madv, 1e-3);     % floor MAD for dims with low variance
Xr = bsxfun(@rdivide, bsxfun(@minus, X0, med), madv); % median-centered, mad-scaled

% 3) PCA reduction
[coeff, score, ~, ~, explained] = pca(Xr);
cumexp = cumsum(explained);
% choose ncomp: first that explains 85-95% (try 90%)
ncomp = find(cumexp > 99.9, 1, 'first');
if isempty(ncomp), ncomp = min(30, size(score,2)); end
Z = score(:, 1:ncomp);

% 3b) Map back to cell per cell structure
idx = cumsum(cellfun(@(x) size(x,1), E));  % row indices
start_idx = [1; idx(1:end-1)+1];
for i = 1:Ncells
    Ti = size(E{i},1);
    Z_cell{i} = Z(start_idx(i):start_idx(i)+Ti-1, :);  % Ti x ncomp
end

% 4) k-means on reduced space
K = 6; % classes
opts = statset('MaxIter',500);
[idx, Cpc] = kmeans(Z, K, 'Replicates', 10, 'Options', opts);

% 4b) Visualize PCA embedding and centroids
figure;
scatter3(Z(:,1), Z(:,2), Z(:,3), 8, idx, 'filled'); 
hold on;
scatter3(Cpc(:,1), Cpc(:,2), Cpc(:,3), 200, 'k', 'filled', 'MarkerEdgeColor','w');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('K-means clusters in PCA space');
grid on; axis equal;

% 5) compute mu_init in original (kept) space
mu_init_kept = zeros(K, sum(keep));
for k = 1:K
    mu_init_kept(k,:) = mean(X0(idx==k,:), 1); % take the mean of all observations with same cluster identity
end

% 6) expand back to full dims, filling dropped dims with global median
mu_init = repmat(median(allE,1), K, 1);
mu_init(:, keep) = mu_init_kept;

% 7) robust s_j for all dims
s_j = mad(allE, 1); % median absolute deviation of each dimension
s_j = max(s_j, 0.04); % floor deviations to avoid normalizing with very small numbers later on

% Shrunken centroid parameters on reduced PCA space (first ~33 comp)
s_j = mad(Z,1) * 10;
s_j = max(s_j,0.02);
mu_init = Cpc;
E = Z_cell';

% Initialize random parameters for beta
% Beta contains the weights for the features influencing the transitions
beta = cell(K,1);
for j = 1:K
    beta{j} = 0.01 * randn(p, K);
end

loglik_prev = -inf;
loglik_trace=[];

%%
% EM LOOP
for iter = 1:maxIter
    fprintf('EM iteration %d\n', iter);

    sum_gamma = zeros(1, K);
    sum_weighted_E = zeros(K, ncomp); % changed from full set of dimensions d to pca comp's ncomp
    total_loglik = 0;
    gamma_all = cell(Ncells,1);
    xi_all = cell(Ncells,1);

    % E-STEP
    % For each cell
    for i = 1:Ncells
        Ti = size(E{i},1); % Lenght of cell track
        logB = zeros(Ti, K);

        % Build emission probability matrix
        % (First from initialized params, then from updated params)
        % For each state/ class K
        for k = 1:K
            % Calculate distance between embedding and shrunken centroid,
            % normalize for comparison
            diffs = (E{i} - mu_init(k,:)) ./ (s_j); % subtract centroid of cluster, normalize by each dims med abs dev
            d2 = sum(diffs.^2, 2); % Squared sum of distances for each timepoint to class K
            logB(:,k) = -0.5 * d2; % log-likelohood of belonging to state at each timepoint
        end

         % Handle missing observations
        missing_observations = propagated_labels{i};
        if(sum(missing_observations)) > 0
            logB(logical(missing_observations)',:) = zeros(sum(missing_observations),K);
        end
        hard_label_missing = hard_labels{i} == 3;
        if(sum(hard_label_missing) > 0)
            logB(hard_label_missing',:) = zeros(sum(hard_label_missing),K); % uniform log-probability across states
        end

        % --- Log-sum-exp normalization for each row ---
        mx = max(logB, [], 2);              % 1×Ti row-wise max
        logB = bsxfun(@minus, logB, mx);    % subtract rowmax from each element in row to prevent overflow
        B = exp(logB); % from log-space to actual (unnormalized) probabilities
        B = bsxfun(@rdivide, B, sum(B,2) + eps); % normalize each row

        % Build transition probability matrix A for each timepoint
        % (First from initialized beta, then from updated beta)
        A = cell(Ti-1,1);

        % For each timepoint, calculate weights of features influencing
        % transitions
        for t = 1:(Ti-1)
            A_t = zeros(K, K);

            % For each state
            for j = 1:K
                logits = X{i}(t,:) * beta{j}; % scores for next states k
                logits = logits - max(logits); 
                Pjk = exp(logits); % to probability
                Pjk = Pjk / sum(Pjk); % normalize
                A_t(j,:) = Pjk; % transition probabilities to states from state j
            end
            A{t} = A_t;
        end
        
        T = size(B,1); % no of timepoints
        K = size(B,2); % no of states

        alpha = zeros(T, K);
        beta_fb = zeros(T, K);
        scale = zeros(T,1);

        % Forward-backward 
        % (remember that this happens for each observed sequence individually)
        
        % Forward pass
        % Do first timepoint separately
        alpha(1, :) = pi_k .* B(1, :); % prob of starting in each state and emitting 1st obs
        scale(1) = sum(alpha(1, :));
        alpha(1, :) = alpha(1, :) / (scale(1) + eps); % normalize, sums to 1

        % For each timepoint
        for t = 2:T
            alpha(t, :) = (alpha(t-1, :) * A{t-1}) .* B(t, :); % the probability of the last timepoint, the transition and the emission
            scale(t) = sum(alpha(t, :));
            alpha(t, :) = alpha(t, :) / (scale(t) + eps);
        end

        % Backward pass
        beta_fb(T, :) = ones(1, K) / (scale(T) + eps); % scale backward probabilities at last timepoint
        % Loop backwards in time
        for t = (T-1):-1:1
            beta_fb(t, :) = (beta_fb(t+1, :) .* B(t+1, :)) * (A{t})'; % probability of observing future data given current state
            beta_fb(t, :) = beta_fb(t, :) / (scale(t) + eps);
        end
        
        loglik_cell_i = sum(log(scale + eps)); % compute log-likelihood of observed sequence

        % Summarize
        total_loglik = total_loglik + loglik_cell_i;

        gamma = alpha .* beta_fb; % Calculate gamma for each state and timepoint
        gamma = gamma ./ sum(gamma, 2);

        xi = cell(Ti-1,1);

        for t = 1:(Ti-1)
            xi_num = (alpha(t,:)' .* A{t}) .* (B(t+1,:) .* beta_fb(t+1,:)); % joint probability of being in state i at time t and state j at time t+1
            %        weighted transit prob    prob of observing future given B(t+1) = k
            xi_num = xi_num / sum(xi_num(:)); % normalize
            xi{t} = xi_num;
        end

        gamma_all{i} = gamma;
        xi_all{i} = xi;

        % For each state k
        for k = 1:K
            w = gamma(:,k); % prob of being in state k at each time point for cell i
            sum_gamma(k) = sum_gamma(k) + sum(w); % summed posterior probability of hidden state equals K for this cell (over whole sequence)
            sum_weighted_E(k, :) = sum_weighted_E(k, :) + (w' * E{i}); % every embedding is weighted
        end

    end

    'Arrived at M-step'

    % M-STEP
    pi_k = zeros(1, K);
    % Loop through the observations
    for i = 1:Ncells
        pi_k = pi_k + gamma_all{i}(1, :); % posterior probability of each state at time t=1 for cell i
                                          % accumulate expected number of sequences starting in each state across all cells
    end
    pi_k = pi_k / sum(pi_k); % normalize
    
    % For each state
    for j = 1:K
        Xtrain = [];
        Wtrain = [];
        
        % For each cell
        for i = 1:Ncells
            Ti = size(E{i},1); % length of sequence

            % Loop over all consecutive pairs of time points
            for t = 1:(Ti-1)
                w_j = sum(xi_all{i}{t}(j,:)); % total probability that the system was in state j at t
                if w_j > 0
                    Xtrain = [Xtrain; X{i}(t,:)]; % Append the features at time t to the training set
                    Wtrain = [Wtrain; xi_all{i}{t}(j,:)]; % append the expected number of transitions
                end
            end
        end
        beta{j} = fit_weighted_multinom(Xtrain, Wtrain, lambda);
        j
    end

    mu_raw = zeros(K, ncomp);
    for k = 1:K
        mu_raw(k, :) = sum_weighted_E(k, :) / (sum_gamma(k) + eps);
    end

    total_weight = sum(sum_gamma);
    mu_global = (sum(sum_weighted_E, 1) / (total_weight + eps));
    delta = (mu_raw - mu_global) ./ (s_j + 1e-8);
    delta_shrunk = sign(delta) .* max(abs(delta) - Delta, 0);
    mu_shrunk = mu_global + delta_shrunk .* s_j; % Uppdaterar vi inte s_j?
    
    % Check for convergence
    if iter > 1 && abs(total_loglik - loglik_prev) < tol
        fprintf('Converged at iteration %d\n', iter);
        break;
    end

    fprintf('Iter %2d | loglik = %.4f | Δloglik = %.4f\n', iter, total_loglik, total_loglik - loglik_prev);
    fprintf('Iter %2d | mean(mu_shrunk)=%.3f | std(s_j)=%.3f | norm(beta)=%.3f\n', ...
        iter, mean(mu_shrunk(:)), std(s_j(:)), norm(cell2mat(beta(:))));
    loglik_trace(iter) = total_loglik;

    loglik_prev = total_loglik;
end

%% SUPPORT FUNCTIONS
function BETA = fit_weighted_multinom(X, W, lambda)
% Fits a multinomial logistic regression model with L2 regularization using
% an iteratively reweighted least squares algorithm

    [N, p] = size(X); % number of samples, predictors
    K = size(W, 2); % number of classes
    BETA = zeros(p, K);
    maxIRLS = 20;

    % For each state
    for k = 1:K
        y = W(:, k);
        beta_k = zeros(p, 1);
        
        for it = 1:maxIRLS
            eta = X * beta_k;
            mu = 1 ./ (1 + exp(-eta));
            S = mu .* (1 - mu);
            z = eta + (y - mu) ./ max(S, 1e-6);
            Wmat = diag(S);
            beta_k = (X' * Wmat * X + lambda * eye(p)) \ (X' * Wmat * z);
        end
        BETA(:, k) = beta_k;
    end
end
