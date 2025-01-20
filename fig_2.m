% Fig 2
% This script produces all the material seen in figure 2 of the paper.
%
% @authors: Madeleine Skeppås
% @date: 15012025
%% Plot confusion matrix 6 class
load('trainedNetwork_6class_v2.mat');

imds = imageDatastore('/Users/madsk418/Desktop/training_data_6class_v2/' ...
       ,"FileExtensions",'.tif', 'LabelSource','foldernames','IncludeSubfolders',true);

[imdsTrain,imdsValidation,imdsTest] = splitEachLabel(imds,0.8,0.1);

imdsTest.ReadFcn = @customReadFunction5;

classNames = categories(imdsTest.Labels);

X = readall(imdsTest);
probs = [];

for i=1:length(X)
    probs = [probs; trainedNetwork_6class_v2.predict(X{i})];
    i
end    

predictedLabels = onehotdecode(probs,classNames,2);
trueLabels = imdsTest.Labels;

confmat = confusionmat(trueLabels,predictedLabels);

cm = confusionchart(confmat, classNames);
cm.Title = 'Classify morphology on test set';
fontsize('scale',1.5)

accuracy = sum(diag(confmat)) / sum(confmat(:));

%% Apply HMM model to tracks
 
% Load table
tbl = load("slice_statistics_03-Sep-2024.mat");
tbl = tbl.slice_statistics;

% Force control samples to have zero dose
tbl.dose(tbl.perturbation == "control") = 0;

% Remove ROIs without tracked objects
idx = cellfun(@(x) isequal(x, []), tbl.traxs);
tbl = tbl(~idx,:);

pool = gcp(); % Start parallel pool

% Fit and apply a Hidden Markov Model to cell tracks
tbl_fit = fit_apply_hmm_v3(tbl);

%% Generate 3D plot, proportions, cell speed, TADs, switches and clustergram

% Fig 2G
plotCellTrajectory(tbl_fit.Viterbi_paths2{1885,1},tbl_fit.traxs{1885,1}, tbl_fit.trays{1885,1})

% Fig 2H
mode = "mean";
morph = "morph";
perturbation = {"PIK 75 HCl"};
varnames = {'Branching','Diffuse translocation','Junk','Locomotion',['Perivascular' ' translocation'],'Round'};
tbl_ext = plot_proportions(tbl_fit, mode, perturbation, varnames,morph);

% Fig 2I, K, L
celline = "U3013MG";
perturbation = "control";
plot_angles_and_speed(tbl_ext, celline, perturbation)

% Fig J
u=unique(tbl_fit.HGCC);
F=[];
for i=1:length(u)
    f1=strmatch('control',tbl_fit.perturbation);
    f2=strmatch(u{i},tbl_fit.HGCC);
    f=intersect(f1,f2);
    F(i,1)=f(1);
end

M=[];

for i=1:length(F)
    M(:,:,i)=tbl_fit.A_est{F(i)};
end

for i=1:length(F)
    subplot(2,3,i);
    plot_transition_hierarchy(min(0.15,M(:,:,i)),varnames,"Junk",14)
    title(u{i})
end
set(gcf,'Position',[100 100 1000 500])
fontsize('scale',1.5)

% Fig M
I=[1 2 4 5 6];
M=M(I,I,:);
M2={};
for i=1:6
    M2{i} = M(:,:,i);
end
cluster_matrices_euclidian(M2, u)



