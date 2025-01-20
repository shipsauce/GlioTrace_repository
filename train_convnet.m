%% create imds
% imds = imageDatastore('/Volumes/MyGroups$/Iron/konfokalmikroskop/Hitesh Montage and Overlays/training_data_man+auto/' ...
%        ,"FileExtensions",'.tif', 'LabelSource','foldernames','IncludeSubfolders',true);

imds = imageDatastore('/Users/madsk418/Desktop/training_data_6class_v2/' ...
       ,"FileExtensions",'.tif', 'LabelSource','foldernames','IncludeSubfolders',true);

imds.ReadFcn = @customReadFunction;

[imdsTrain,imdsValidation,imdsTest] = splitEachLabel(imds,0.8,0.1);

%% perform zero mean centering

train_images = readall(imdsTrain);

image = zeros(size(train_images{1}));

for i=1:length(train_images)
    im = train_images{i};
    image = image + im;
    i
end

mean_image = image / length(train_images);
% save('mean_image_6class.mat', 'mean_image');

inputSize = size(readimage(imds,1));

%% perform oversampling of training set
labels=imdsTrain.Labels;
[G,classes] = findgroups(labels);
numObservations = splitapply(@numel,labels,G);

desiredNumObservationsPerClass = max(numObservations);

files = splitapply(@(x){randReplicateFiles(x,desiredNumObservationsPerClass)},imdsTrain.Files,G);
files = vertcat(files{:});
labels=[];
info=strfind(files,'/');
for i=1:numel(files)
    idx=info{i};
    dirName=files{i};
    targetStr=dirName(idx(end-1)+1:idx(end)-1);
    targetStr2=cellstr(targetStr);
    labels=[labels;categorical(targetStr2)];
end
imdsTrain.Files = files;
imdsTrain.Labels=labels;
labelCount_oversampled = countEachLabel(imdsTrain);

%% adjust customReadFunction to include zero mean centering operation, perform augmentation

imdsTrain.ReadFcn = @customReadFunction2;
imdsTest.ReadFcn = @customReadFunction2;
imdsValidation.ReadFcn = @customReadFunction2;

augmenter = imageDataAugmenter( ...
RandXReflection=true, ...
RandRotation=[-90 90], ...
RandScale=[0.5 2]);

augimdsTrain = augmentedImageDatastore(inputSize(1:2),imdsTrain,DataAugmentation=augmenter);
augimdsValidation = augmentedImageDatastore(inputSize(1:2),imdsValidation, DataAugmentation=augmenter);

numClasses = 5;
%% construct convnet
layers = [
    imageInputLayer(inputSize, 'Name', 'input', 'Normalization', 'none')

    convolution2dLayer(3, 32, 'Padding', 'same', 'Name', 'conv1')
    reluLayer('Name', 'relu1')
    batchNormalizationLayer('Name', 'bn1')

    maxPooling2dLayer(2, 'Stride', 2, 'Name', 'pool1')

    convolution2dLayer(3, 64, 'Padding', 'same', 'Name', 'conv2')
    reluLayer('Name', 'relu2')
    batchNormalizationLayer('Name', 'bn2')

    maxPooling2dLayer(2, 'Stride', 2, 'Name', 'pool2')

    convolution2dLayer(3, 128, 'Padding', 'same', 'Name', 'conv3')
    reluLayer('Name', 'relu3')
    batchNormalizationLayer('Name', 'bn3')

    maxPooling2dLayer(2, 'Stride', 2, 'Name', 'pool3')

    convolution2dLayer(3, 256, 'Padding', 'same', 'Name', 'conv4')
    reluLayer('Name', 'relu4')
    batchNormalizationLayer('Name', 'bn4')

    fullyConnectedLayer(512, 'Name', 'fc1')
    reluLayer('Name', 'relu5')

    fullyConnectedLayer(256, 'Name', 'fc2')
    reluLayer('Name', 'relu6')

    fullyConnectedLayer(numClasses, 'Name', 'output')
    softmaxLayer('Name', 'softmax')
];

%% sanity check with mnist dataset
% dataFolder = fullfile(toolboxdir('nnet'),'nndemos','nndatasets','DigitDataset');
% imds = imageDatastore(dataFolder, ...
%     'IncludeSubfolders',true,'LabelSource','foldernames');
% 
% numTrainFiles = 750;
% [imdsTrain,imdsValidation] = splitEachLabel(imds,numTrainFiles,"randomize");
% 
% img = readimage(imds,1);
% inputSize = size(img);
% 
% numClasses = length(unique(imds.Labels));
% 
% options = trainingOptions("sgdm", ...
%     InitialLearnRate=0.01, ...
%     MaxEpochs=4, ...
%     Shuffle="every-epoch", ...
%     ValidationData=imdsValidation, ...
%     ValidationFrequency=30, ...
%     Plots="training-progress", ...
%     Metrics="accuracy", ...
%     Verbose=false);
% 
% net = trainnet(imdsTrain,layers,"crossentropy",options);

%% sanity check with small training set

numImagesToUse = 20;

subsetIndices = randperm(length(imdsTrain.Files), numImagesToUse);
subsetIndices_valid = randperm(length(imdsValidation.Files), numImagesToUse);

imdsTrainSubset = subset(imdsTrain, subsetIndices);
imdsValidationSubset = subset(imdsValidation, subsetIndices_valid);

augmenter = imageDataAugmenter( ...
RandXReflection=true, ...
RandRotation=[-90 90], ...
RandScale=[0.5 2]);

train_sub = augmentedImageDatastore(inputSize(1:2),imdsTrainSubset,DataAugmentation=augmenter);
valid_sub = augmentedImageDatastore(inputSize(1:2),imdsValidationSubset, DataAugmentation=augmenter);

miniBatchSize = 5;

numObservationsTrain = numel(imdsTrainSubset.Files);
numIterationsPerEpoch = floor(numObservationsTrain/miniBatchSize);

options = trainingOptions('adam', ...
    MiniBatchSize=miniBatchSize, ...
    MaxEpochs=50, ...
    InitialLearnRate=0.001, ...
    LearnRateSchedule='piecewise', ...
    ValidationData=valid_sub, ...
    OutputNetwork="best-validation-loss", ...
    Plots="training-progress", ...
    ValidationFrequency=5, ...
    Verbose=false, ...
    L2Regularization=0, ...
    Shuffle="every-epoch", ...
    Metrics="accuracy");


net = trainnet(train_sub,layers, 'crossentropy',options);

%% load net
load('lgraph_1.mat')
net = dlnetwork(lgraph_1);
plot(net)

%% define training process

options = trainingOptions("sgdm", ...
    InitialLearnRate=0.01, ...
    MaxEpochs=30, ...
    Shuffle="every-epoch", ...
    ValidationData=augimdsValidation, ...
    ValidationFrequency=30, ...
    Plots="training-progress", ...
    Metrics="accuracy", ...
    Verbose=true, ...
    MiniBatchSize=64, ...
    LearnRateSchedule="piecewise", ...
    L2Regularization=1e-2);

%% train network
net_trained = trainnet(augimdsTrain,net,"crossentropy",options);

%% evaluate

classNames = categories(imdsTest.Labels);

X = readall(imdsTest);
probs = [];

for i=1:length(X)
    probs = [probs; net_trained.predict(X{i})];
    i
end    

predictedLabels = onehotdecode(probs,classNames,2);
trueLabels = imdsTest.Labels;

confmat = confusionmat(trueLabels,predictedLabels);

cm = confusionchart(confmat, classNames);
cm.Title = 'Classify morphology on test set';
fontsize('scale',1.5)

accuracy = sum(diag(confmat)) / sum(confmat(:))