function A_log = evaluate_models_on_trajectories(models, A_raw, trajectory, K, eps)
% Returns (T, K, K) array of log transition probabilities.
%
%
% @authors: André Lasses Armatowski, Linnea Hallin 
% @date: 09112025

if nargin < 4
    eps = 1e-6;
end

[T, F] = size(trajectory); % Timepoints x features of single observation
A_log = zeros(T, K, K);

if isempty(models)
    % Use raw counts, nornalize convert to log space
    A_raw_norm = A_raw ./ (sum(A_raw,2) + eps);
    A_log_base = log(A_raw_norm + eps);

    % Repeat for all timesteps
    A_log = repmat(A_log_base, 1, 1, T);
    A_log = permute(A_log, [3 1 2]);   % (T, K, K)
else
    for i = 1:numel(models)
        % trajectory: (T, F)
        % model.coef: (K, F)
        % model.intercept: (K,)
        %logits = trajectory * models{i}.coef.' + repmat(models{i}.intercept(:).', T, 1); % compute model weighted features
        [~,logits] = models{i}.predict(trajectory);
        A_log(:, i, :) = log_softmax_mat(logits, 2); % convert into probabilities in log space
    end
end

end
