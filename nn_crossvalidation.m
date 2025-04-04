% 10 fold cv
% split into train test validation 10x --> titta på spridning av accuracy
% data = 80% train, 10% test, 10% valid

imds = imageDatastore('/Users/madsk418/Desktop/training_data_6class_v2/' ...
       ,"FileExtensions",'.tif', 'LabelSource','foldernames','IncludeSubfolders',true);

imds.ReadFcn = @customReadFunction;

load('lgraph_6.mat')

metric = [];

inputSize = size(readimage(imds,1));

% Stratified split into 10 folds
[fold1, fold2, fold3, fold4, fold5, fold6, fold7, fold8, fold9, fold10] = splitEachLabel(imds,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1);
            
folds = {fold1.Files, fold2.Files, fold3.Files, fold4.Files, fold5.Files, fold6.Files, fold7.Files, fold8.Files, fold9.Files, fold10.Files};

fold_indices = [1:8, 9, 10;
                2:9, 10, 1;
                3:10, 1, 2;
                [4:10 1], 2, 3;
                [5:10 1:2], 3, 4;
                [6:10 1:3], 4, 5;
                [7:10 1:4], 5, 6;
                [8:10 1:5], 6, 7;
                [9:10 1:6], 7, 8;
                [10 1:7], 8, 9];

% Perform cross validation of network performance
for i=1:10

    % Training data
    training_folds = fold_indices(i, 1:8);
    files_training = [];
    for j=1:8
        files_training = [files_training; folds{training_folds(j)}];
    end

    imdsTrain = imageDatastore(files_training, "IncludeSubfolders",true, "LabelSource","foldernames");

    % Validation data
    imdsValidation = imageDatastore(folds{fold_indices(i, 9)}, "IncludeSubfolders",true, "LabelSource","foldernames");

    % Test data
    imdsTest = imageDatastore(folds{fold_indices(i, 10)}, "IncludeSubfolders",true, "LabelSource","foldernames");

    % Perform oversampling of training set
    labels=imdsTrain.Labels;
    [G,classes] = findgroups(labels);
    numObservations = splitapply(@numel,labels,G);
    
    desiredNumObservationsPerClass = max(numObservations);
    
    files = splitapply(@(x){randReplicateFiles(x,desiredNumObservationsPerClass)},imdsTrain.Files,G);
    files = vertcat(files{:});
    labels=[];
    info=strfind(files,'/');
    for p=1:numel(files)
        idx=info{p};
        dirName=files{p};
        targetStr=dirName(idx(end-1)+1:idx(end)-1);
        targetStr2=cellstr(targetStr);
        labels=[labels;categorical(targetStr2)];
    end

    imdsTrain.Files = files;
    imdsTrain.Labels=labels;
    labelCount_oversampled = countEachLabel(imdsTrain);

    % Adjust customReadFunction to include zero mean centering operation, perform augmentation

    imdsTrain.ReadFcn = @customReadFunction5;
    imdsTest.ReadFcn = @customReadFunction5;
    imdsValidation.ReadFcn = @customReadFunction5;
    
    augmenter = imageDataAugmenter( ...
    RandXReflection=true, ...
    RandRotation=[-90 90], ...
    RandScale=[0.5 2]);
    
    augimdsTrain = augmentedImageDatastore(inputSize(1:2),imdsTrain,DataAugmentation=augmenter);
    augimdsValidation = augmentedImageDatastore(inputSize(1:2),imdsValidation, DataAugmentation=augmenter);
    
    numClasses = 6;
    
    net = dlnetwork(lgraph_6);

    % Define training process

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

    % Train network
    net_trained = trainnet(augimdsTrain,net,"crossentropy",options);

    classNames = categories(imdsTest.Labels);

    % Evaluate on test set
    X = readall(imdsTest);
    probs = [];
    
    for m=1:length(X)
        probs = [probs; net_trained.predict(X{m})];
        m
    end    
    
    predictedLabels = onehotdecode(probs,classNames,2);
    trueLabels = imdsTest.Labels;
    
    confmat = confusionmat(trueLabels,predictedLabels);
    
    accuracy = sum(diag(confmat)) / sum(confmat(:));

    metric(i) = accuracy;

end

% metric_full = [0.8784    0.8939    0.9091    0.9080    0.8849    0.8825    0.8947    0.9259    0.9399    0.9130];

%% TME NN 10-fold cv-validation

imds = imageDatastore('/Volumes/MyGroups$/Iron/konfokalmikroskop/Hitesh Montage and Overlays/training_data_tme/v2/' ...
       ,"FileExtensions",'.tif', 'LabelSource','foldernames','IncludeSubfolders',true);

imds.ReadFcn = @customReadFunction;

load('lgraph_2class.mat')

metric2 = [];

inputSize = size(readimage(imds,1));

% Stratified split into 10 folds
[fold1, fold2, fold3, fold4, fold5, fold6, fold7, fold8, fold9, fold10] = splitEachLabel(imds,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1);
            
folds = {fold1.Files, fold2.Files, fold3.Files, fold4.Files, fold5.Files, fold6.Files, fold7.Files, fold8.Files, fold9.Files, fold10.Files};

fold_indices = [1:8, 9, 10;
                2:9, 10, 1;
                3:10, 1, 2;
                [4:10 1], 2, 3;
                [5:10 1:2], 3, 4;
                [6:10 1:3], 4, 5;
                [7:10 1:4], 5, 6;
                [8:10 1:5], 6, 7;
                [9:10 1:6], 7, 8;
                [10 1:7], 8, 9];

% Perform cross validation of network performance
for i=1:10

    % Training data
    training_folds = fold_indices(i, 1:8);
    files_training = [];
    for j=1:8
        files_training = [files_training; folds{training_folds(j)}];
    end

    imdsTrain = imageDatastore(files_training, "IncludeSubfolders",true, "LabelSource","foldernames");

    % Validation data
    imdsValidation = imageDatastore(folds{fold_indices(i, 9)}, "IncludeSubfolders",true, "LabelSource","foldernames");

    % Test data
    imdsTest = imageDatastore(folds{fold_indices(i, 10)}, "IncludeSubfolders",true, "LabelSource","foldernames");

    % Perform oversampling of training set
    labels=imdsTrain.Labels;
    [G,classes] = findgroups(labels);
    numObservations = splitapply(@numel,labels,G);
    
    desiredNumObservationsPerClass = max(numObservations);
    
    files = splitapply(@(x){randReplicateFiles(x,desiredNumObservationsPerClass)},imdsTrain.Files,G);
    files = vertcat(files{:});
    labels=[];
    info=strfind(files,'/');
    for p=1:numel(files)
        idx=info{p};
        dirName=files{p};
        targetStr=dirName(idx(end-1)+1:idx(end)-1);
        targetStr2=cellstr(targetStr);
        labels=[labels;categorical(targetStr2)];
    end

    imdsTrain.Files = files;
    imdsTrain.Labels=labels;
    labelCount_oversampled = countEachLabel(imdsTrain);

    % Adjust customReadFunction to include zero mean centering operation, perform augmentation

    imdsTrain.ReadFcn = @customReadFunction3;
    imdsTest.ReadFcn = @customReadFunction3;
    imdsValidation.ReadFcn = @customReadFunction3;
    
    augmenter = imageDataAugmenter( ...
    RandXReflection=true, ...
    RandRotation=[-90 90], ...
    RandScale=[0.5 2]);
    
    augimdsTrain = augmentedImageDatastore(inputSize(1:2),imdsTrain,DataAugmentation=augmenter);
    augimdsValidation = augmentedImageDatastore(inputSize(1:2),imdsValidation, DataAugmentation=augmenter);
    
    numClasses = 2;
    
    net = dlnetwork(layers_1);

    % Define training process

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

    % Train network
    net_trained = trainnet(augimdsTrain,net,"crossentropy",options);

    classNames = categories(imdsTest.Labels);

    % Evaluate on test set
    X = readall(imdsTest);
    probs = [];
    
    for m=1:length(X)
        probs = [probs; net_trained.predict(X{m})];
        m
    end    
    
    predictedLabels = onehotdecode(probs,classNames,2);
    trueLabels = imdsTest.Labels;
    
    confmat = confusionmat(trueLabels,predictedLabels);
    
    accuracy = sum(diag(confmat)) / sum(confmat(:));

    metric2(i) = accuracy;

end

% metric2 = [0.9362    0.9833    0.9630    0.9732    0.9832    0.9532    0.9564    0.9866    0.9262    0.9765];
