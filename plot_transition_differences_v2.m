function plot_transition_differences_v2(before_matrices, after_matrices, state_names, plot_titles)
% This function takes two state transition matrices and plots a directed
% graph which is the difference between the first (control) and second 
% (treatment) matrix.
%
% Input parameters:
% before_matrices, after_matrices: cell arrays of NxN matrices representing
% transition rates before and after treatment, respectively.
% state_names: A cell array of strings containing custom names for the states.
% plot_titles: A cell array of strings containing titles for each subplot.
%
% @date: 20012025    


% before_matrices, after_matrices: cell arrays of NxN matrices representing
% transition rates before and after treatment, respectively.
% state_names: A cell array of strings containing custom names for the states.
% plot_titles: A cell array of strings containing titles for each subplot.

    % Validate input cell arrays
    num_matrices = numel(before_matrices);
    if num_matrices ~= numel(after_matrices)
        error('Number of before-treatment and after-treatment matrices must be the same');
    end
    
    % Determine the number of states and validate each matrix
    [num_states, ~] = size(before_matrices{1});
    for i = 1:num_matrices
        if any(size(before_matrices{i}) ~= [num_states, num_states]) || ...
           any(size(after_matrices{i}) ~= [num_states, num_states])
            error('All matrices must be square and of the same dimensions');
        end
    end
    
    % Validate state names
    if length(state_names) ~= num_states
        error('The number of state names must match the number of states in the matrices');
    end
    
    % Remove the third state (both from matrices and names)
    for i = 1:num_matrices
        before_matrices{i}(3, :) = [];
        before_matrices{i}(:, 3) = [];
        after_matrices{i}(3, :) = [];
        after_matrices{i}(:, 3) = [];
    end
    state_names(3) = []; % Remove the third state name
    num_states = num_states - 1; % Update state count after removal
    
    % Compute all difference matrices
    diff_matrices = cell(1, num_matrices);
    max_weight = 0; % To track the max weight across all difference matrices
    for i = 1:num_matrices
        diff_matrices{i} = after_matrices{i} - before_matrices{i};
        diff_matrices{i}(1:num_states+1:end) = 0; % Set diagonal to zero
        max_weight = max(max_weight, max(abs(diff_matrices{i}(:))));
    end
    
    % Determine layout for subplots (e.g., 2x2, 3x2, etc.)
    num_cols = ceil(sqrt(num_matrices));
    num_rows = ceil(num_matrices / num_cols);

    % Plot each difference matrix as a subplot with consistent scaling and color
    for i = 1:num_matrices
        G = digraph(diff_matrices{i});
        
        % Edge properties
        edge_weights = G.Edges.Weight;
        
        % Find the top 10 strongest edges by absolute weight
        [~, sorted_idx] = sort(abs(edge_weights), 'descend');
        top10_idx = sorted_idx(1:min(10, length(sorted_idx))); % Indices of top 10 strongest edges
        
        % Initialize edge colors to white and set thickness based on top 10 edges only
        edge_colors = repmat([1, 1, 1], height(G.Edges), 1); % Default color is white
        edge_thickness = ones(height(G.Edges), 1); % Default thickness is minimal
        
        % Determine color and thickness for each of the top 10 edges
        max_thickness = 10; % Scaling for edge thickness
        for j = 1:length(top10_idx)
            idx = top10_idx(j);
            weight = edge_weights(idx);
            if weight > 0
                % Upregulated (red)
                edge_colors(idx, :) = [1, 0, 0] * min(abs(weight) / max_weight, 1);
                edge_colors(idx, :) = [1, 0, 0];
            else
                % Downregulated (blue)
                edge_colors(idx, :) = [0, 0, 1] * min(abs(weight) / max_weight, 1);
            end
            % Scale thickness
            edge_thickness(idx) = 1 + (max_thickness - 1) * (abs(weight) / max_weight);
        end
        
        % Plot directed graph with customized properties
        h = plot(G, 'Layout', 'layered', 'EdgeColor', edge_colors, 'ArrowSize', 15, 'NodeFontSize', 15, 'NodeLabel', state_names);
        h.NodeColor = [0.7, 0.7, 0.7]; % Nodes are gray
        h.LineWidth = edge_thickness;
        h.MarkerSize = 25;
        title(plot_titles)
        
    end
    
end
