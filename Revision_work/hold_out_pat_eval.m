%% Hold-out model training to evaluate patient-specific performance
% Load data
imds = imageDatastore('/Users/madsk418/Desktop/training_data_6class_v2/' ...
       ,"FileExtensions",'.tif', 'LabelSource','foldernames','IncludeSubfolders',true);

imds.ReadFcn = @customReadFunction;

load('lgraph_6.mat')

metric = [];

inputSize = size(readimage(imds,1));

metadata=readtable('/Volumes/MyGroups$/Iron/konfokalmikroskop/Hitesh Montage and Overlays/hitesh_metadata.xlsx');

%% Fetch patient celline labels
patlabs = categorical([]);

for i=1:length(imds.Files)
    % Convert to numeric if needed
    try
        expid = str2double(regexp(imds.Files{i}, 'exp(\d{3})', 'tokens', 'once'));
        try
            hgcc = metadata.HGCC{expid};
        catch
            expid = str2double(regexp(imds.Files{i}, 'exp(\d{2})', 'tokens', 'once'));
            hgcc = metadata.HGCC{expid};
        end
        hgcc = categorical(cellstr(hgcc(1:7)));
        patlabs(i) = hgcc;
    catch
        if(contains(imds.Files{i}, "score"))
            expr = '(?<=stack_)\d+';
            stackid = regexp(imds.Files{i}, expr, 'match');
            expid = stack_to_exp_mapping_autogen_snapshots_v2(str2num(stackid{:}));
            hgcc = metadata.HGCC{expid};
            hgcc = categorical(cellstr(hgcc(1:7)));
            patlabs(i) = hgcc;
        else
            expr = '(?<=stack_)\d+';
            stackid = regexp(imds.Files{i}, expr, 'match');
            try
                expid = stack_to_exp_mapping_autogen_snapshots_v1(str2num(stackid{:}));
                hgcc = metadata.HGCC{expid};
                hgcc = categorical(cellstr(hgcc(1:7)));
                patlabs(i) = hgcc;
            catch
                patlabs(i) = categorical(NaN);
            end
        end
    end
    i
end

pat = unique(patlabs');
pat = pat(1:5);
pat_count = length(pat);

%% Hold-out loop: in each loop, hold out another patient celline and evaluate performance after training

for i=1:pat_count
    idx = (patlabs == pat(i))';

    imdsHoldout = subset(imds, ~idx);
    imdsHeldout = subset(imds, idx);

    [imdsTrain,imdsValidation] = splitEachLabel(imdsHoldout, 0.8, 0.2);

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
    imdsHeldout.ReadFcn = @customReadFunction5;
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

    classNames = categories(imdsHeldout.Labels);

    % Evaluate on test set
    X = readall(imdsHeldout);
    probs = [];
    
    for m=1:length(X)
        probs = [probs; net_trained.predict(X{m})];
        m
    end    
    
    predictedLabels = onehotdecode(probs,classNames,2);
    trueLabels = imdsHeldout.Labels;
    
    confmat = confusionmat(trueLabels,predictedLabels);
    
    accuracy = sum(diag(confmat)) / sum(confmat(:));

    metric{i} = confmat;

end

% metric =
% 
%     0.7077    0.8305    0.9108    0.9714    0.9103

%% Per class recall
recall = {};
idx = [1 2 3 5];
% Calculate per class recall for each patient
% For each patient
for i=1:length(pat)-1
    confmat = metric{idx(i)};
    
    % For each class
    for j=1:6
        calc = confmat(j,j)/ (sum(confmat(j,:)));
        try
            recall{j} = [recall{j} calc];
        catch
            recall{j} = calc;
        end
    end
end

% Plot distribution as boxplot
colors = [
     24 120 187;  % Electric Blue
    239  71 111;  % Hot Pink Red
    255 196  37;  % Lemon Yellow
     38 201 164;  % Aqua Green
    155  89 182;  % Vivid Purple
    255 126  52   % Bright Orange
] / 255;

mat = cell2mat(recall');
figure;
for i=1:4
    scatter(1:6, mat(:,i), 300, colors(i,:), "filled","diamond", "DisplayName", string(pat(idx(i))))
    hold on
end
ax=gca;
xticks(1:6)
xticklabels({"Branching", "DT", "Crowded", "Locomotion", "PT", "Round"})
ax.XLim = [0 7];
ax.YLim = [0 1.1];
ax.Box = 1;
ylabel('Per class recall')
title("Patient-wise per class recall")
fontsize('scale',1.6)
legend("Location","southeast")

%% Färgfemman

% Patient-normalized mean conf mat
matagg = zeros(6,6);
max2mean = [];
matraw = [];
idx = [1 2 3 5];
for i=1:4
    matnorm = metric{idx(i)} ./ max(metric{idx(i)}, [], 2);
    matnorm(isnan(matnorm)) = 0;
    matagg = matagg + matnorm;
    max2mean(:,:,i) = matnorm;
    matraw(:,:,i) = metric{idx(i)};
end
matagg = matagg./4;

max2mean_rat = max(max2mean, [], 3) ./ nanmean(max2mean, 3);

classNames = {"Branching", "DT", "Crowded", "Locomotion", "PT", "Round"};

imagesc(matagg)
colormap cool
title('Patient-aggregated confusion matrix')
xticks(1:6)
yticks(1:6)
xticklabels(classNames)
yticklabels(classNames)
fontsize('scale', 1.6)




