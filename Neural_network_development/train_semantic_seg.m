% This script defines training data and a U-net architecture for semantic segmentation. 
% 5-fold cross validation is performed followed by evaluation on a held out
% test set.
%
% @authors: Madeleine Skeppås
% @date: 20022025

% Load data and hold out test set
dataDir = '/Users/madsk418/UU Dropbox/Madeleine S/Simulation_and_invasion/comp/output/Madeleine/Semantic_seg_quarters/ims_manualseg_only/';
maskDir = '/Users/madsk418/UU Dropbox/Madeleine S/Simulation_and_invasion/comp/output/Madeleine/Semantic_seg_quarters/labs_manualseg_only/';

imds = imageDatastore(fullfile(dataDir), 'FileExtensions', '.png');
classes = [0 1];
labs = pixelLabelDatastore(fullfile(maskDir), ["background" "vasc"], classes);

% Split data using subset
trainRatio = 0.9;
testRatio =  0.1;
numImages = length(imds.Files);
indices = randperm(numImages);
trainInd = indices(1:round(trainRatio*numImages));
testInd = indices(round((trainRatio)*numImages)+1:end);

imdsTrain = subset(imds, trainInd);
imdsTest = subset(imds, testInd);

pxdsTrain = subset(labs, trainInd);
pxdsTest = subset(labs, testInd);

tbl = countEachLabel(labs);
imageFreq = tbl.PixelCount ./ tbl.ImagePixelCount;
classWeights = median(imageFreq) ./ imageFreq;

%%

% Define U-Net layers
inputSize = [250 250 1];  % Image size
numClasses = 2;           % Binary segmentation (background & vessels)
numFilters = 32;          % Base number of filters (kept low for small dataset)

% Input Layer
layers = [
    imageInputLayer(inputSize, 'Normalization', 'none', 'Name', 'input')
    
    % Encoder (Downsampling Path)
    convolution2dLayer(3, numFilters, 'Padding', 'same', 'Name', 'conv1_1')
    reluLayer('Name', 'relu1_1')
    convolution2dLayer(3, numFilters, 'Padding', 'same', 'Name', 'conv1_2')
    reluLayer('Name', 'relu1_2')
    maxPooling2dLayer(2, 'Stride', 2, 'Name', 'pool1')
    
    convolution2dLayer(3, numFilters*2, 'Padding', 'same', 'Name', 'conv2_1')
    reluLayer('Name', 'relu2_1')
    convolution2dLayer(3, numFilters*2, 'Padding', 'same', 'Name', 'conv2_2')
    reluLayer('Name', 'relu2_2')
    maxPooling2dLayer(2, 'Stride', 2, 'Name', 'pool2')
    
    convolution2dLayer(3, numFilters*4, 'Padding', 'same', 'Name', 'conv3_1')
    reluLayer('Name', 'relu3_1')
    convolution2dLayer(3, numFilters*4, 'Padding', 'same', 'Name', 'conv3_2')
    reluLayer('Name', 'relu3_2')
    maxPooling2dLayer(2, 'Stride', 2, 'Name', 'pool3')
    
    convolution2dLayer(3, numFilters*8, 'Padding', 'same', 'Name', 'conv4_1')
    reluLayer('Name', 'relu4_1')
    dropoutLayer(0.3, 'Name', 'drop4') % Dropout to prevent overfitting
    convolution2dLayer(3, numFilters*8, 'Padding', 'same', 'Name', 'conv4_2')
    reluLayer('Name', 'relu4_2')
    
    % Decoder (Upsampling Path)
    transposedConv2dLayer(2, numFilters*4, 'Stride', 2, 'Name', 'upsample3')
    convolution2dLayer(3, numFilters*4, 'Padding', 'same', 'Name', 'conv5_1')
    reluLayer('Name', 'relu5_1')
    convolution2dLayer(3, numFilters*4, 'Padding', 'same', 'Name', 'conv5_2')
    reluLayer('Name', 'relu5_2')
    
    transposedConv2dLayer(2, numFilters*2, 'Stride', 2, 'Name', 'upsample2')
    convolution2dLayer(3, numFilters*2, 'Padding', 'same', 'Name', 'conv6_1')
    reluLayer('Name', 'relu6_1')
    convolution2dLayer(3, numFilters*2, 'Padding', 'same', 'Name', 'conv6_2')
    reluLayer('Name', 'relu6_2')
    
    transposedConv2dLayer(2, numFilters, 'Stride', 2, 'Name', 'upsample1')
    convolution2dLayer(3, numFilters, 'Padding', 'same', 'Name', 'conv7_1')
    reluLayer('Name', 'relu7_1')
    convolution2dLayer(3, numFilters, 'Padding', 'same', 'Name', 'conv7_2')
    reluLayer('Name', 'relu7_2')
    
    % Final Convolution Layer
    convolution2dLayer(1, numClasses, 'Padding', 'same', 'Name', 'conv_final')
    softmaxLayer('Name', 'softmax')
];

net = dlnetwork(layers);

%% CV model

options = trainingOptions('sgdm', ...
    'MaxEpochs', 120, ...
    'MiniBatchSize', 12, ...
    'Shuffle', 'every-epoch', ...
    'ValidationFrequency', 5, ...
    'Verbose', true, ...
    'OutputNetwork','best-validation-loss', ...
    'Metrics','fscore', ...
    'ValidationPatience',30,...         
    'InitialLearnRate',0.001);

f1scores = cv_eval(net, options, imdsTrain, pxdsTrain, classWeights);

%% Extract scores
validation_scores=[];
for i=1:5
    [best_loss,row] = min(f1scores{i,2}.Loss);
    validation_scores(i) = f1scores{i,2}.FScore(row);
end

%% Assess performance on test set

numTest = numel(imdsTest.Files);
f1Scores = zeros(numTest, 1);

for i = 1:numTest
    % Read Image & Label
    testImg = readimage(imdsTest, i);
    testLabel = readimage(pxdsTest, i);
    
    % Predict segmentation mask using semanticseg
    predictedMask = semanticseg(testImg, net_trained);
    predictedMask = ~(predictedMask == 'C1');

    groundTruthMask = ~(testLabel == 'background');
    
    % Display results
    subplot(2,2,1); imshow(testImg); title('Test Image');
    subplot(2,2,2); imshow(predictedMask); title('Predicted Segmentation');
    subplot(2,2,3); imshowpair(testImg, predictedMask)
    subplot(2,2,4); imshow(~(testLabel == 'background'))
    hold on

    tp = sum(predictedMask(:) & groundTruthMask(:)); % True Positives
    fp = sum(predictedMask(:) & ~groundTruthMask(:)); % False Positives
    fn = sum(~predictedMask(:) & groundTruthMask(:)); % False Negatives
    
    precision = tp / (tp + fp + eps);
    recall = tp / (tp + fn + eps);
    f1Scores(i) = 2 * (precision * recall) / (precision + recall + eps);
end
