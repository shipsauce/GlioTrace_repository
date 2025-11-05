%% Shrunken centroid init based on labelled data

load("trainedNetwork_6class_v2.mat")

imds = imageDatastore('/Users/madsk418/Desktop/training_data_6class_v2/' ...
       ,"FileExtensions",'.tif', 'LabelSource','foldernames','IncludeSubfolders',true);

imds.ReadFcn = @customReadFunction;

classNames = categories(imds.Labels);

X = readall(imds);
probs = [];
embeddings_relu6 = [];
embeddings_fc2 = [];

% Get softmax probs and 256 embedding for each image in labelled data

for i=1:length(X)
    probs = [probs; trainedNetwork_6class_v2.predict(X{i})];
    feature_embedding = predict(trainedNetwork_6class_v2,X{i},'Outputs', 'relu6');
    embeddings_relu6 = [embeddings_relu6; feature_embedding];
    embeddings_fc2 = [embeddings_fc2; predict(trainedNetwork_6class_v2,X{i},'Outputs', 'fc2')];
    i
end

trueLabels = imds.Labels;

classes = unique(trueLabels);

dims = size(embeddings_relu6,2);

centroids = zeros(6,dims);

for i=1:length(classes)
    idx = trueLabels == classes(i);
    centroids(i,:) = mean(embeddings_relu6(idx,:),1);
end

dim_std = std(embeddings_relu6,1);

% Sanity check

[coeff, score, ~, ~, explained] = pca(embeddings_relu6);

centroids_pca = centroids * coeff(:, 1:3);  % Project onto the first 2 principal components

figure;
scatter3(score(:,1), score(:,2), score(:,3), 8,trueLabels,'filled'); % Plot embeddings in PCA space
hold on;
scatter3(centroids_pca(:,1), centroids_pca(:,2), centroids_pca(:,3) ,100, 'filled', 'r'); % Plot centroids
title('PCA Projection of Embeddings and Centroids');
xlabel('PC1');
ylabel('PC2');

% Compare spread of centroids to spread of data

var_data = var(embeddings_relu6, 0, 1);  % Variance per dimension across all samples
std_centroids = std(centroids, 0, 1);  % Standard deviation across centroids

figure;
plot(var_data, 'b', 'LineWidth', 2);
hold on;
plot(std_centroids, 'r--', 'LineWidth', 2);
legend('Data Variance', 'Centroid Std Dev');
xlabel('Dimension');
ylabel('Value');

% Save

% Concatenate centroids and standard deviations into one matrix
output_data = [centroids; dim_std];

% Save to text file
writematrix(output_data, 'shrunken_centroids_init.txt', 'Delimiter', ' ');






























