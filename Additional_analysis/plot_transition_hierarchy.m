function plot_transition_hierarchy(transition_matrix, state_names, exclude_state, show_number)
% This function takes a state transition matrix and plots a directed
% graph.
%
% Input parameters:
% transition_matrix - A NxN matrix representing transition rates between states.
% state_names - A cell array of strings containing custom names for the states.
% exclude_state - A string specifying the state name to exclude from the network.
% show_number - The number of edges to display in the plot
%
% @date: 16012025

% Validate the dimensions of the transition matrix and state names
[num_states, cols] = size(transition_matrix);
if num_states ~= cols
    error('The transition matrix must be square (NxN)');
end

if length(state_names) ~= num_states
    error('The number of state names must match the number of states in the transition matrix');
end

% Find the index of the state to exclude
exclude_idx = find(strcmp(state_names, exclude_state));
if isempty(exclude_idx)
    error('The specified state to exclude does not exist in state_names');
end

% Remove the row and column corresponding to the excluded state from the transition matrix
transition_matrix(exclude_idx, :) = [];
transition_matrix(:, exclude_idx) = [];

% Remove the corresponding state name
state_names(exclude_idx) = [];

% Exclude self-loops by setting diagonal elements to zero
num_remaining_states = size(transition_matrix, 1);
transition_matrix(1:num_remaining_states+1:end) = 0; % Set diagonal elements to zero

% Create a directed graph object from the modified transition matrix
G = digraph(transition_matrix);

colors = [0.8660,    0.3290,  0;
0.3290,    0.7130,    1.0000;
0.9960,    0.5640,    0.2620;
0.4540,    0.9210,    0.8540;
     0,   0.6390,   0.6390
];

% Calculate edge thickness based on transition rates with dramatic scaling
min_thickness = 1;
max_thickness = 100*max(transition_matrix(:)); % Maximum thickness for the edges
edge_weights = G.Edges.Weight;
scaled_weights = (edge_weights / max(edge_weights)).^2; % Square the normalized weights for dramatic scaling
edge_thickness = min_thickness + (max_thickness - min_thickness) * scaled_weights;

% Find the 10 strongest transitions
[~, sorted_idx] = sort(edge_weights, 'descend');
top10_idx = sorted_idx(1:min(show_number, length(sorted_idx)));

% Create a color vector for edges
edge_colors = repmat([1 1 1], height(G.Edges), 1); % Default color is white
edge_colors(top10_idx, :) = repmat([0.5 0.5 0.5], length(top10_idx), 1); % Top 10 edges in grey

% Create a label vector for edges
edge_labels = repmat({''}, height(G.Edges), 1); % Default label is empty
weights = round(G.Edges.Weight(top10_idx)./min(G.Edges.Weight(top10_idx)),2);
edge_labels(top10_idx) = cellstr(num2str(weights)); % Set labels for top 10 edges

% Plot the transition network with the top 10 strongest transitions and individual colors
h = plot(G, 'Layout', 'layered', 'EdgeColor', edge_colors, 'ArrowSize', 15,'EdgeFontSize',35,'NodeFontSize',15);

% Customize appearance of the plot
 h.NodeLabel = [];  % Use remaining state names for the nodes
h.NodeColor = colors; % Assign individual colors to nodes
h.LineWidth = edge_thickness; % Set edge thickness based on dramatic scaling of transition rates
h.MarkerSize = 35; % Increase the size of the node markers
ax = gca;
end
