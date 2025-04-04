function clusters = cluster_matrices_euclidian(matrices, labels)
% This function takes a cell array of state transition matrices and cluster
% them based on the Euclidian distance of the non-diagonal elements using
% hierarchical clustering.
% 
% Input parameters:
% matrices - A cell array where each cell contains a matrix.
% labels - A cell array of labels for the matrices (for x-axis ticks).
%
% @date: 16012025

% Number of matrices
numMatrices = length(matrices);

% Ensure the number of labels matches the number of matrices
if nargin < 2 || length(labels) ~= numMatrices
    error('Number of labels must match the number of matrices.');
end

% Initialize distance matrix
distMatrix = zeros(numMatrices);

% Compute pairwise Euclidean distance between matrices after removing diagonals
for i = 1:numMatrices
    for j = i+1:numMatrices
        % Remove diagonal from matrices
        matrix1 = matrices{i};
        matrix2 = matrices{j};
        matrix1(logical(eye(size(matrix1)))) = NaN; % Set diagonal elements to NaN
        matrix2(logical(eye(size(matrix2)))) = NaN; % Set diagonal elements to NaN

        % Flatten matrices into vectors, excluding NaN (diagonal) values
        vec1 = matrix1(:);
        vec2 = matrix2(:);
        validIdx = ~isnan(vec1) & ~isnan(vec2); % Only use non-NaN elements

        % Compute Euclidean distance between the two vectors
        euclideanDist = norm(vec1(validIdx) - vec2(validIdx));

        % Store the Euclidean distance in the distance matrix
        distMatrix(i, j) = euclideanDist;
        distMatrix(j, i) = euclideanDist;  % Symmetric matrix
    end
end

% Perform hierarchical clustering using the distance matrix
Z = linkage(squareform(distMatrix), 'complete');

% Automatically determine clusters using an inconsistency threshold
inconsistencyMatrix = inconsistent(Z);
threshold = mean(inconsistencyMatrix(:, 4)) + std(inconsistencyMatrix(:, 4));

% Form clusters using the computed threshold
clusters = cluster(Z, 'cutoff', threshold, 'criterion', 'distance');

% Plot dendrogram with the cutoff line and colored clusters
figure;
[h, T, perm] = dendrogram(Z, 0, 'ColorThreshold', threshold); % Color by the threshold

for i=1:5
    h(i).LineWidth = 1.5;
    h(i).Color = [0 0 0];
end

ax=gca;
ax.Box = 1;

% Set the x-axis labels using the input labels
xticks(1:numMatrices); % Ensure the x-ticks correspond to each matrix
xticklabels(labels(perm)); % Permute the labels to match the dendrogram order
xtickangle(45); % Rotate labels for better readability if necessary

% Plot the threshold line on the dendrogram
hold on;
line([0 numMatrices+1], [threshold threshold], 'Color', 'r', 'LineStyle', '--'); % Plot the threshold line
title('Hierarchical Clustering Dendrogram with Auto Cutoff (Euclidean Distance)');
xlabel('Matrix Labels');
ylabel('Euclidean Distance');

% Display the clusters for debugging or information purposes (optional)
disp('Cluster assignments:');
disp(clusters);
end

