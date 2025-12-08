function models = update_transitions(trajectories, xis, glm_iters, penalty, prev_models)
% trajectories: cell array of (T, F)
% xis:          cell array of (T, K, K)
% Returns:      cell array of K multinomial logistic models
% @authors: André Lasses Armatowski, Linnea Hallin 
% @date: 09112025

if nargin < 3
    glm_iters = 500;
end
if nargin < 4
    penalty = 1;
end
if nargin < 5
    prev_models = [];
end

% Determine K from first xi
first_xi = xis{1};
K = size(first_xi, 2);

models = cell(1, K);

total_obs = sum(cellfun(@(x) (size(x,1)-1)*K, xis));

for i = 1:K
    X = zeros(total_obs, size(trajectories{1},2));
    y = zeros(total_obs,1);
    sample_weights = zeros(total_obs,1);
    idx = 1;
    N_seq = numel(xis);
    
    for c = 1:N_seq
        xi = xis{c};
        feats = trajectories{c};
        T = size(xi, 1);

        tt = 1:(T-1);
        idx2 = numel(tt)*K;
        
        X(idx:idx+idx2-1,:) = repelem(feats(tt,:), K, 1);
        y(idx:idx+idx2-1) = repmat((1:K)', numel(tt),1);
        tmp = squeeze(xi(tt,i,1:K));
        sample_weights(idx:idx+idx2-1) = reshape(tmp',[],1);

        % for t = 1:T-1
        %     for j = 1:K
        %         X = [X; feats(t, :)];
        %         y = [y; j];
        %         sample_weights = [sample_weights; xi(t, i, j)];
        %     end
        % end
        idx = idx + idx2;
        
    end

    opts = statset('MaxIter', glm_iters);

    % Test vanilla HMM
    T = array2table(X);
    T.Y = y;
    mdl = fitmnr(T, 'Y ~ 1', 'Weights', sample_weights);

    % mdl = fitmnr(X, y, ...
    % 'Weights', sample_weights);

    models{i} = mdl;
end
end
