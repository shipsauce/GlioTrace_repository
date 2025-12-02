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

for i = 1:K
    X = [];
    y = [];
    sample_weights = [];

    N_seq = numel(xis);
    for c = 1:N_seq
        xi = xis{c};
        feats = trajectories{c};
        T = size(xi, 1);

        for t = 1:T-1
            for j = 1:K
                X = [X; feats(t, :)];
                y = [y; j];
                sample_weights = [sample_weights; xi(t, i, j)];
            end
        end
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
