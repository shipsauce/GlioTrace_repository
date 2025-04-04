function val_f1 = cv_eval(net, options, trainval_data, trainval_labs, classWeights)
% This function takes a neural network architecture, a trainingOptions
% object and data for training and validation, and performs 5-fold cross
% validation of the network.
%
% Input parameters:
% net - neural network architecture
% options - trainingOptions object
% trainval_data - training data as an imageDatastore object
% trainval_labs - training labels as a pixelDatastore object
% classWeights - weights calculated based on the prevalence of pixels from
%       the different classes, used in the crossentropy calculation
%
% Output parameters:
% val_f1 - table of training- and validation history for each (fifth)
%            iteration for each fold
%
% @authors: Madeleine Skeppås
% @date: 17012025


% List files
ims = trainval_data.Files;
labs = trainval_labs.Files;

num_files = length(ims);

val_f1 = [];

% Define partitions
indices = randperm(num_files); % Randomly shuffle indices
fold_size = floor(num_files / 5); % Base size of each fold
remainder = mod(num_files, 5); % Extra files to distribute

fold_indices = cell(5, 2);
start_idx = 1;

for i = 1:5
    % Determine the number of elements in this fold
    extra = i <= remainder; % Distribute remainder among first few folds
    num_in_fold = fold_size + extra;
    
    % Assign indices to fold
    fold_indices{i,2} = indices(start_idx : start_idx + num_in_fold - 1);
    fold_indices{i,1} = setdiff(indices, indices(start_idx : start_idx + num_in_fold - 1));
    
    % Update start index
    start_idx = start_idx + num_in_fold;
end

% CV loop
for i=1:5
    fprintf(['CV Fold: ' num2str(i) '\n'])
    % Create PixelLabelImageDatastores for training and validation
    imds_train = subset(trainval_data,fold_indices{i,1});
    pxds_train = subset(trainval_labs,fold_indices{i,1});

    imds_val = subset(trainval_data,fold_indices{i,2});
    pxds_val = subset(trainval_labs,fold_indices{i,2});

    trainingData = pixelLabelImageDatastore(imds_train,pxds_train);

    validationData = pixelLabelImageDatastore(imds_val,pxds_val);

    options.ValidationData = validationData;

    % Fit weights
    [net_trained, info] = trainnet(trainingData, net, @(Y,T) modelLoss(Y,T,classWeights), options);

    % Extract info on validation accuracy of best model based on best
    % validation loss
    val_f1{i,1} = info.TrainingHistory;
    val_f1{i,2} = info.ValidationHistory;
end

end

%%

function loss = modelLoss(Y,T,classWeights)
    weights = dlarray(classWeights,"C");
    mask = ~isnan(T);
    T(isnan(T)) = 0;
    loss = crossentropy(Y,T,weights,Mask=mask,NormalizationFactor="mask-included");
end