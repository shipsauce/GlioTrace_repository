function [slice_statistics, vasculature_statistics] = build_statistics_v3(stackfile, output)
% 
% Given a set of ROIs (stacks) from a brain slice culture experiment, 
% this script will calculate ROI-level and cell-level statistics and 
% summarize it in a stacktable with information on all ROIs.
% 
% ROI-level statistics: 
% - growth rate 
% - total green signal in every frame
% - simplified SAD (sum of absolute differences) 
% - adMAD (adjusted mean of absolute differences)
%
% Cell-level statistics: 
% - traX: x-coordinates of individual cells across frames
% - traY: y-coordinates of individual cells across frames
% - props: information on cell shape & TME interactions across frames
%
% @authors: Madeleine Skeppås, Sven Nelander
% @date: 05062024
%
% Load neural networks and create stacktable
metadata=readtable('/Volumes/MyGroups$/Iron/konfokalmikroskop/Hitesh Montage and Overlays/hitesh_metadata.xlsx');

load('brainslice_manuscript_repo/trainedNetwork_6class_v2.mat');
load('brainslice_manuscript_repo/trainedNetwork_tme.mat');
blocksize=61;
warning('off', 'all')

stacktable=build_stack_table(metadata,stackfile);

% Define the set of stacks to be analyzed

cellines= unique(stacktable.HGCC); % Retrieve the names of HGCC cellines

perturbations = "all";
if strcmp(perturbations, "all")
    idx = true(height(stacktable), 1); % Select all rows
else
    idx = ismember(stacktable.perturbation, perturbations);
end

% Save the number of input arguments to toggle visualisation later
args = nargin;

subtable = stacktable(idx,:);

% Calculate statistics

% Initialize arrays to store data
growth_rate=[];
sad=[];
adMAD = {};
sum_green = {};
msd_curves = {};
traxs = {};
trays = {};
props = {};
vasc_length_stack = {};
segmented_stack = {};
vasculature_statistics = table;

if(nargin == 1)
    output = [];
end

pool = gcp(); % Start parallel pool

net_const = parallel.pool.Constant(trainedNetwork_6class_v2);
tme_net_const = parallel.pool.Constant(trainedNetwork_tme);

fprintf('Building slice statistics table...\n')
fprintf(['Estimated time: ' num2str(ceil((height(subtable) * 22.45)/60)) ' minutes\n'])

% Iterate through the stacks
parfor i=1:height(subtable)
    empty_video = false;
    fprintf(['Tracking cells in stack: ' num2str(i) ' / ' num2str(height(subtable)) '...\n']);
    
    % Load the stack
    stack = load(subtable.file{i}); 
    stack = stack.stack;

    % Save green and red channel in separate variables
    gbm=stack.Tstack;
    vasc=stack.Vstack;
    dt=subtable.delta_t(i);

    % Calculate ROI-level statistics
    
    growth_rate(i,1)=stackscore_naive_growth_rate(gbm,dt);  % Calculate growth rate
    sad(i,1)=stackscore_naive_SAD(gbm); % Calculate simplified SAD
    
    ad_stack = abs(gbm(:,:,2:end)-gbm(:,:,1:end-1)); % Calculate absolute differences
    c = zeros(size(gbm,3)-1,1);
    
    for j=1:size(ad_stack,3) % Iterate over absolute differences
        ad =  ad_stack(:,:,j);
        ad(ad<20) = 0;  % Threshold the image
        c(j) = mean(mean(ad))/mean(mean(mean(gbm(:,:,j:j+1)))); % Calculate adMAD 
    end

    adMAD{i} = c; % Save to variable

    sum_green{i} = squeeze(sum(sum(gbm))); % Calculate the total green signal in every frame

    % Calculate cell-level statistics

    % Set parameters for cell detection
    sigmah = 16; % Average cell size
    hsizeh = 60; % Size of LoG filter used in blob detection
    cutoff = 2e-4; % Intensity threshold in blob detection

    % Identify cell coordinates
    try
        [cellsx,cellsy,intensity,feat,vascc]=macro_track2(gbm,vasc,sigmah,hsizeh,cutoff, blocksize, "normal");
    catch
        try
            % If cells are sparse in frames, try another set of parameters
            [cellsx,cellsy,intensity,feat,vascc]=macro_track2(gbm,vasc,sigmah/2,hsizeh/2,cutoff*6, blocksize, "sparse");
            fprintf(['Sparse ROI no: ' num2str(i) '\n'])
        catch
            empty_video = true;
            fprintf(['Empty ROI no: ' num2str(i) '\n'])
        end
    end
    
    if(~empty_video)
        % Classify cell morphology and interactions with the TME
        properties = classify_tumor_cells(feat, vascc, blocksize, net_const.Value, tme_net_const.Value, i);
    

         % Track cells with a Kalman filter
        [traX, traY, x_hat_history, phenotypes, phenonames, Kn_history, startidx] = track_tumor_cells2(cellsx, cellsy, properties);
    
        % Connect tracklets of fragmented tracks
        [traX, traY, phenotypes] = connect_tracklets(traX, traY, phenotypes);
  
    else
        traX = [];
        traY = [];
        phenotypes = [];

        empty_video = false;
    end
    
    % Save to variables
    traxs{i} = traX;
    trays{i} = traY;
    props{i} = phenotypes;


    [vasc_length, segstack] = segment_quantify_vasculature(vasc,subtable(i,:), output);
    vasc_length_stack(i) = {vasc_length};
    segmented_stack(i) = {segstack};

    % Visualize tracking and classification if output path is given
    if(args > 1)
        mode = "morphology";
        path = output;
        vis_tracking_HD(traX,traY,gbm, vasc,subtable(i,[1 2 3 8 9]), phenotypes, mode, startidx,path);
    
        mode = "tme";
        vis_tracking_HD(traX,traY,gbm, vasc,subtable(i,[1 2 3 8 9]), phenotypes, mode, startidx,path);
    end
end

% Merge stacktable with calculated statistics
vasculature_statistics.vasc_length_stack = vasc_length_stack';
vasculature_statistics.segmented_stack = segmented_stack';
vasculature_statistics = [subtable vasculature_statistics];

subtable.traxs = traxs';
subtable.trays = trays';
subtable.props = props';
subtable.growth_rate = growth_rate;
subtable.sad = sad;
subtable.adMAD = adMAD';
subtable.sum_green = sum_green';

slice_statistics = subtable;
end

