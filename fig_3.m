%% Fig 3
% This script produces material seen in figure 3 of the paper.
%
% @authors: Madeleine Skeppås
% @date: 15012025
%% Quantify vasculature
% Fetch metadata about experiments
metadata=readtable('/Volumes/MyGroups$/Iron/konfokalmikroskop/Hitesh Montage and Overlays/hitesh_metadata.xlsx');

% Build table of ROIs from file
stacktable=build_stack_table(metadata,'vasc_stacks_complete.txt');
subtable = stacktable(stacktable.perturbation == "control",:);

% Exclude video with corrupted intensity in beginning
subtable(31,:) = [];

subtable = sortrows(subtable,"HGCC","ascend");
stacktable_const = parallel.pool.Constant(subtable);

% Create empty variables for storing time-resolved and average
% vessel measurements
vessel_length = {};
vessel_length_average = [];

% Start parallel pool
pool = gcp(); 

%%
% Loop through the image stacks
load("net_2025_02_20.mat")

for i=1:height(subtable)
    % Load the stack
    stack = load(stacktable_const.Value.file{i});
    stack = stack.stack;

    % Save green and red channel in separate variables
    gbm=stack.Tstack;
    vasc=stack.Vstack;

    info = stacktable_const.Value(i,:);
    
    % Quantify vasculature length across stack (red channel)
    [vasc_length_stack,~] = quantify_vasculature_v3(vasc,info, net_trained);

    % Save to variables
    vessel_length{i} = vasc_length_stack;
    vessel_length_average(i) = mean(cell2mat(vasc_length_stack));
    i
end

% Append variables to table
subtable.vessel_length = vessel_length';
subtable.vessel_length_average = vessel_length_average';

%% Save table
vessel_statistics = subtable;
save(['vessel_statistics_' char(datetime('today')) '.mat'], 'vessel_statistics')

%% Fig 3B 
% Segment tif stacks for figure

folderpath='/Users/madsk418/UU Dropbox/Madeleine S/Simulation_and_invasion/comp/output/Madeleine/fig_3_seg';
fps = 10;
video_size_factor=10; % increase for higher resolution video, default = 1

% Specify folder containing TIFF stacks
folder = folderpath;

% Get a list of all TIFF files in the folder
fileList = dir(fullfile(folder, '*.tif'));

% Specify video filename and properties
frames_per_second = fps; % Adjust frame rate as needed

tif_stack = cell(length(fileList),1);
tumor = [];
vasc = [];

% Loop through each TIFF stack
for i = 1:length(fileList)

    % Create VideoWriter object
    movie_filename = [folder '/segstack_2_' fileList(i).name];
    myvideo = VideoWriter(movie_filename, 'MPEG-4');
    myvideo.FrameRate = frames_per_second;
    open(myvideo);

    % Read the TIFF stack
    filename = fullfile(folder, fileList(i).name);
    tif_info = imfinfo(filename); % Get TIFF file information
    num_frames = numel(tif_info); % Get the number of frames in the stack
    tif_single = imread(filename);
    tif = zeros([size(tif_single) num_frames]);

    % Loop through the frames and save channels separately
    for t=1:num_frames
         current_frame = imread(filename, t, 'Info', tif_info);
         tif(:,:,:,t) = current_frame;
         tumor(:,:,t) = current_frame(:,:,2);
         vasc(:,:,t) = current_frame(:,:,1);
    end

    info = 'tif_stack';

    % Segment vasculature
    sensitivity = 0.5;
    [~, segstack] = quantify_vasculature_v3(vasc,info,net_trained);

    % Save segmented images also as a tif-stack
    tifstack_filename = [folder '/segstack_tif_' fileList(i).name];
    im = segstack(:,:,:,1);
    im = uint8(imresize(im, 'scale', video_size_factor));
    imwrite(im, tifstack_filename, 'WriteMode', 'overwrite', 'Compression', 'none');

    for j=2:num_frames
        im = segstack(:,:,:,j);
        imshow(im/255);

        % Convert current frame to video frame
        current_frame = im2frame(uint8(imresize(im, 'scale', video_size_factor)));
    
        % Write current frame to video
        writeVideo(myvideo, current_frame);
        
        % Write to tif-stack
        im = uint8(imresize(im, 'scale', video_size_factor));
        imwrite(im, tifstack_filename, 'WriteMode', 'append', 'Compression', 'none');
    end

    % Close video writer and all figures
    close(myvideo);
    close all force;
end

%% Fig 3C 
% Make time-resolved plot of the relative change (values are set-wise normalized to control)

% Load table
load('vessel_statistics_semseg.mat')
subtable = vessel_statistics;

cellines = unique(subtable.HGCC); % Retrieve the names of HGCC cellines

% Loop through the cell lines
for i=1:length(cellines)
    
    subplot(3,2,i)
    hold on
  
    hgcc = cellines{i};

    % Tab is a subset of subtable only containing stacks from the same
    % cell line
    tab = subtable(strcmp(subtable.HGCC,hgcc),:);

    % Save stacks from tumor vs control region in separate variables (used
    % for calculating 
    tumor_tab = tab(logical(tab.vessel_origin == "tumor"),:);
    ctrl_tab = tab(logical(tab.vessel_origin == "control"),:);
    
    % Find all unqiue sets in this cell line
    sets = unique(tab.set);

    % Create empty variables to store data from all stacks
    values_ctrl = {};
    values_tumor = {};
    
    % Loop over the sets
    for k=1:length(unique(tab.set))
        sett = sets(k);

        % Set_tab only contains stacks from the same set
        set_tab = tab(tab.set == sett,:);

        % Separate set stacks based on vessel status
        tumor = set_tab(logical(set_tab.vessel_origin == "tumor"),:);
        ctrl = set_tab(logical(set_tab.vessel_origin == "control"),:);
        
        % Pad cell arrays with NaNs to achieve equal length of arrays
        ctrl_out = cellfun(@(x) cell2mat(x),ctrl.vessel_length,'UniformOutput',false);
        max_length = max(cellfun(@length, ctrl_out));
        ctrl_filled = cellfun(@(x) [x, NaN(1, max_length - length(x))], ctrl_out, 'UniformOutput', false);

        % Calculate control mean for set
        ctrl_mean = boxFilter(nanmean(cell2mat(ctrl_filled),1),8);

        % Pad cell arrays with NaNs to achieve equal length of arrays
        tumor_out = cellfun(@(x) cell2mat(x),tumor.vessel_length,'UniformOutput',false);
        max_length = max(cellfun(@length, tumor_out));
        tumor_filled = cellfun(@(x) [x, NaN(1, max_length - length(x))], tumor_out, 'UniformOutput', false);

        % Set-wise normalization of values achieved by dividing with
        % control mean
        tumor_norm = cellfun(@(x) x./ctrl_mean, tumor_filled, 'UniformOutput',false);
        ctrl_norm = cellfun(@(x) x./ctrl_mean, ctrl_filled, 'UniformOutput',false);
        
        % Append set curves to variable for entire celline
        values_ctrl = [values_ctrl; ctrl_norm];
        values_tumor = [values_tumor; tumor_norm];
        
    end

    % Create empty variables for storing values
    ctrl_change = nan(height(values_ctrl),max(cellfun(@length, values_ctrl)));
    tumor_change = nan(height(values_tumor),max(cellfun(@length, values_tumor)));
    
    % Loop through values
    for m=1:height(values_ctrl)
        vals = values_ctrl{m}'; % Extract values from cell array
        vals = vals./median(vals(1:5)) * 100; % Convert to procentual change
        ctrl_change(m,1:length(vals)) = vals;
    end
    
    % Do the same for curves from tumor region
    for m=1:height(values_tumor)
        vals = values_tumor{m}'; 
        vals = vals./median(vals(1:5)) * 100; 
        tumor_change(m,1:length(vals)) = vals;
    end
    
    % Plot
    values_ctrl_normalized = ones(size(tumor_change))*100;
    p1 = stdshade_v1(ctrl_change, 0.1, [0.5 0.7 1], ['Control (n = ' num2str(length(unique(ctrl_tab.exp))) ')' ' (ROIs = ' num2str(height(values_ctrl)) ')'],8);
    p2 = stdshade_v1(tumor_change, 0.1, [1 0.251 0.251], ['Tumor (n = ' num2str(length(unique(tumor_tab.exp))) ')' ' (ROIs = ' num2str(height(values_tumor)) ')'],8);
    legend([p1 p2],"Location","best")
    
    % Add figure details
    xlabel('t')
    ylabel('Vessel length (% change)')
    title(hgcc)
end

% Adjust axes
lims = [67 55  60 45 65 100];
for i=1:6
    subplot(3,2,i)
    a=axis();
    a(2)= lims(i);
    axis(a)
    ax = gca;
    ax.YAxis.Limits = [40 120];
    ax.Box=1;
end

% Add subplot title and scale
sgtitle('Vasculature length in area near tumor vs normal vasculature')
fontsize('scale',1.5)

%% Fig 3D 
% Make violin plot of AUC of time-resolved change based on set-wise normalized curves (tumor vs control)

% Load table
load('vessel_statistics_semseg.mat')
subtable = vessel_statistics;

% Retrieve the names of HGCC cellines
cellines = unique(subtable.HGCC);

tumor_auc = nan(78,6);
control_auc = nan(78,6);
group = [];
data={};
c = [0.5 0.7 1;
    1 0.251 0.251];

% Loop through the cell lines
for i=1:length(cellines)
    hgcc = cellines{i};

    % Tab is a subset of subtable only containing stacks from the same
    % cell line
    tab = subtable(strcmp(subtable.HGCC,hgcc),:);

    % Save stacks from tumor vs control region in separate variables (used
    % for calculating 
    tumor_tab = tab(logical(tab.vessel_origin == "tumor"),:);
    ctrl_tab = tab(logical(tab.vessel_origin == "control"),:);
    
    % Find all unqiue sets in this cell line
    sets = unique(tab.set);

    % Create empty variables to store data from all stacks
    values_ctrl = {};
    values_tumor = {};
    
    % Loop over the sets
    for k=1:length(unique(tab.set))
        sett = sets(k);

        % Set_tab only contains stacks from the same set
        set_tab = tab(tab.set == sett,:);

        % Separate set stacks based on vessel status
        tumor = set_tab(logical(set_tab.vessel_origin == "tumor"),:);
        ctrl = set_tab(logical(set_tab.vessel_origin == "control"),:);
        
        % Pad cell arrays with NaNs to achieve equal length of arrays
        ctrl_out = cellfun(@(x) cell2mat(x),ctrl.vessel_length,'UniformOutput',false);
        max_length = max(cellfun(@length, ctrl_out));
        ctrl_filled = cellfun(@(x) [x, NaN(1, max_length - length(x))], ctrl_out, 'UniformOutput', false);

        % Calculate control mean for set
        ctrl_mean = boxFilter(nanmean(cell2mat(ctrl_filled),1),8);

        % Pad cell arrays with NaNs to achieve equal length of arrays
        tumor_out = cellfun(@(x) cell2mat(x),tumor.vessel_length,'UniformOutput',false);
        max_length = max(cellfun(@length, tumor_out));
        tumor_filled = cellfun(@(x) [x, NaN(1, max_length - length(x))], tumor_out, 'UniformOutput', false);

        % Set-wise normalization of values achieved by dividing with
        % control mean
        tumor_norm = cellfun(@(x) x./ctrl_mean, tumor_filled, 'UniformOutput',false);
        ctrl_norm = cellfun(@(x) x./ctrl_mean, ctrl_filled, 'UniformOutput',false);
        
        % % Append set curves to variable for entire celline
        values_ctrl = [values_ctrl; ctrl_norm];
        values_tumor = [values_tumor; tumor_norm];
        
    end

    % Create empty variables for storing values
    ctrl_change = nan(height(values_ctrl),max(cellfun(@length, values_ctrl)));
    tumor_change = nan(height(values_tumor),max(cellfun(@length, values_tumor)));
    
    % Loop through values
    for m=1:height(values_ctrl)
        vals = values_ctrl{m}'; % Extract values from cell array
        vals = vals./median(vals(1:5)) * 100; % Convert to procentual change
        ctrl_change(m,1:length(vals)) = vals;
    end
    
    % Do the same for curves from tumor region
    for m=1:height(values_tumor)
        vals = values_tumor{m}'; 
        vals = vals./median(vals(1:5)) * 100;
        tumor_change(m,1:length(vals)) = vals;
    end
    
   % Create empty variables for storing AUC values
   tc_auc_celline = [];
   ctrl_auc_celline = [];

   % Calculate AUC for each curve and normalize by length
   for p=1:size(tumor_change,1)
        tc_auc_celline = [tc_auc_celline trapz(tumor_change(p,~isnan(tumor_change(p,:))))/sum(~isnan(tumor_change(p,:)))];

        ctrl_auc_celline = [ctrl_auc_celline trapz(ctrl_change(p,~isnan(ctrl_change(p,:))))/sum(~isnan(ctrl_change(p,:)))];
   end

   % Merge with data from all other cellines
   tumor_auc(1:length(tc_auc_celline),i) = tc_auc_celline;
   control_auc(1:length(ctrl_auc_celline),i) = ctrl_auc_celline;
end

% Reformat for violin plot
data{2} = tumor_auc;
data{1} = control_auc;

% Plot
h = daviolinplot(data,'colors', c,'outsymbol','k+',...
    'boxcolors','same','scatter',1,'jitter',1, 'scattercolors', 'same','xtlabels', cellines,...
    'legend',{'Control', 'Tumor'}, 'scatteralpha', 0.5, 'violinalpha', 0.7);
fontsize('scale',1.5)

ax=gca;
ax.Box = 1;
ax.YLim = [-10 220];

%% Calculate significance from AUC values
p_values=[];
for i=1:length(cellines)
    p = anova1([data{1}(:,i) data{2}(:,i)]);
    p_values = [p_values p];
end

%%  Fig 3F, G
% Calculate vasculature association and microglia association statistics

% Load table
tbl = load("slice_statistics_03-Sep-2024.mat");
tbl = tbl.slice_statistics;

% Force control samples to have zero dose
tbl.dose(tbl.perturbation == "control") = 0;

% Remove ROIs without tracked objects
idx = cellfun(@(x) isequal(x, []), tbl.traxs);
tbl = tbl(~idx,:);

% Plot TME label proportions
style = "mean";
perturbation = {"control"};
label = "tme";
varnames = {'Microglia colocalized', 'Vessel associated', 'Non-associated'};
tbl_ext = plot_proportions(tbl, style, perturbation, varnames, label);

% Select relevant ROIs
tbl_ext_ctrl = tbl_ext(tbl_ext.perturbation == "control",:);

tbl_ext_ctrl.Microglia_interaction_events = cell(height(tbl_ext_ctrl),1);
tbl_ext_ctrl.Microglia_interaction_time = cell(height(tbl_ext_ctrl),1);
tbl_ext_ctrl.Norm_MIE= cell(height(tbl_ext_ctrl),1);

% Iterate through the stacktable and add information on microglia
% interaction events (frequency and duration)
for i=1:height(tbl_ext_ctrl)
    row = tbl_ext_ctrl(i,:);
    tme_stat = row.props{1}{8}; % Get TME label matrix (time x cells)
    
    interaction_events = [];
    interaction_time = [];

    % For each cell...
    for j=1:width(tme_stat)
        path = tme_stat(:,j);
        
        idx = extractSubsequences(path, 1); % ... extract the subsequences where label is 1 (Mg assoc)
        if(~isempty(idx))
            % If the sequence contains at least one such sequence, add +1 to the counter
            interaction_events = [interaction_events 1]; 
        end
        
        % ... count the number of frames that the interaction lasted for each
        % subsequence
        interaction_time = [interaction_time cellfun(@(x) length(x), idx)];
    end

    % Summarize the number of interaction events and the differing
    % interaction lengths for this particular ROI
    tbl_ext_ctrl.Microglia_interaction_events(i) = {sum(interaction_events)};
    tbl_ext_ctrl.Microglia_interaction_time(i) = {(interaction_time-1) .* tbl_ext_ctrl.delta_t(i)};
    tbl_ext_ctrl.Norm_MIE(i) = {(sum(interaction_events)/width(tme_stat))*100}; % Percentage of cells interacting with Mg out of all cells in ROI
    i
end

% Calculate the average interaction time per ROI
tbl_ext_ctrl.Average_microglia_interac_time = cellfun(@(x) mean(x), tbl_ext_ctrl.Microglia_interaction_time, UniformOutput=false);

tbl_ext_ctrl.HGCC = categorical(tbl_ext_ctrl.HGCC); 
tbl_ext_ctrl.Average_microglia_interac_time = cell2mat(tbl_ext_ctrl.Average_microglia_interac_time);
tbl_ext_ctrl.Microglia_interaction_events = cell2mat(tbl_ext_ctrl.Microglia_interaction_events);
tbl_ext_ctrl.Norm_MIE = cell2mat(tbl_ext_ctrl.Norm_MIE);

col = [
  [115/255, 191/255, 184/255];    
  [115/255, 191/255, 184/255];   
  [115/255, 191/255, 184/255];  
  [115/255, 191/255, 184/255];   
  [115/255, 191/255, 184/255];
  [115/255, 191/255, 184/255]
];

% Plot violin plot of percentage of cells interacting with microglia
figure;
daviolinplot(tbl_ext_ctrl.Norm_MIE,'groups', tbl_ext_ctrl.HGCC,'colors', col,'outsymbol','k+',...
    'boxcolors','same','scatter',1,'jitter',1, 'scattercolors', 'same',...
    'scatteralpha', 0.4, 'violinalpha', 1, 'boxalpha', 0.9,'xtlabels', unique(tbl_ext_ctrl.HGCC));
ax=gca;
ax.YLim = [-40 160];
fontsize('scale',1.3)
ax.Box = 1;

% Perform ANOVA
[p, tbl, stats] = anova1(tbl_ext_ctrl.Norm_MIE, tbl_ext_ctrl.HGCC);
multcompare(stats)

% Plot boxplot of average interaction time
figure;
b = boxchart(tbl_ext_ctrl.HGCC, tbl_ext_ctrl.Average_microglia_interac_time);
ax=gca;
ax.YLim = [0 80];
b.SeriesIndex = 7;
ax.Box = 1;
fontsize('scale',1.3)

%% Plot confusion matrix of TME network

load('trainedNetwork_tme.mat')

imds = imageDatastore('/Volumes/MyGroups$/Iron/konfokalmikroskop/Hitesh Montage and Overlays/training_data_tme/v2', ...
    "FileExtensions",'.tif', 'LabelSource','foldernames','IncludeSubfolders',true);

[imdsTrain,imdsValidation,imdsTest] = splitEachLabel(imds,0.8,0.1);

imdsTest.ReadFcn = @customReadFunction3;

classNames = categories(imdsTest.Labels);

X = readall(imdsTest);
probs = [];

for i=1:length(X)
    probs = [probs; trainedNetwork_tme.predict(X{i})];
    i
end    

predictedLabels = onehotdecode(probs,classNames,2);
trueLabels = imdsTest.Labels;

confmat = confusionmat(trueLabels,predictedLabels);

cm = confusionchart(confmat, classNames);
cm.Title = 'Classify TME on test set';
fontsize('scale',1.5)

accuracy = sum(diag(confmat)) / sum(confmat(:));


%%
function dataOut = boxFilter(dataIn, fWidth)
% apply 1-D boxcar filter for smoothing
fWidth = fWidth - 1 + mod(fWidth,2); %make sure filter length is odd
dataStart = cumsum(dataIn(1:fWidth-2),2);
dataStart = dataStart(1:2:end) ./ (1:2:(fWidth-2));
dataEnd = cumsum(dataIn(length(dataIn):-1:length(dataIn)-fWidth+3),2);
dataEnd = dataEnd(end:-2:1) ./ (fWidth-2:-2:1);
dataOut = conv(dataIn,ones(fWidth,1)/fWidth,'full');
dataOut = [dataStart,dataOut(fWidth:end-fWidth+1),dataEnd];
end
