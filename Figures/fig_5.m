%% Fig 5
% This script produces material seen in figure 5 of the paper.
%
% @authors: Madeleine Skeppås
% @date: 15012025
%% Build ROI table of stainings from textfile
roitab=readtable('fig5_stainings_roi.txt','ReadVariableNames',false,'delimiter','\t');
roitab.Properties.VariableNames{1} = 'Folderpath';
roitab.ROI_no = zeros(height(roitab),1);
roitab.HGCC = cell(height(roitab),1);
roitab.Yellow_ch = cell(height(roitab),1);
roitab.Area = cell(height(roitab),1);
roitab.Perturbation = cell(height(roitab),1);
roitab.Dose = zeros(height(roitab),1);

%% Create metadata

invalid_rows = [];

for i=1:height(roitab)

    try
        % ROI number
        roitab.ROI_no(i) = str2double(regexp(roitab.Folderpath{i}, '(\d+)(?=\.tif)', 'match'));
    
        % HGCC
        hgcc = char(regexp(roitab.Folderpath{i}, '(?<=BS\d+-)\d+(?=-)', 'match'));

        if(isempty(hgcc))
            hgcc = char(regexp(roitab.Folderpath{i}, '(?<=BS\d+-)\d+(?=_)', 'match'));
        end
        
        roitab.HGCC{i} = ['U' hgcc 'MG'];
    
        % Yellow channel cell population
        roitab.Yellow_ch{i} = char(regexp(roitab.Folderpath{i}, '(?<=set\d+_)\w+(?=_)', 'match'));
        
        % Core/ edge
        if(contains(roitab.Folderpath{i}, 'edge', 'IgnoreCase',true))
            roitab.Area{i} = 'edge';
        else
            roitab.Area{i} = 'core';
        end
    
        % Perturbation
        if(contains(roitab.Folderpath{i}, 'cont',"IgnoreCase",true))
            roitab.Perturbation{i} = 'control';
        elseif(contains(roitab.Folderpath{i}, 'dasa',"IgnoreCase",true))
            roitab.Perturbation{i} = 'dasatinib';
        elseif(contains(roitab.Folderpath{i}, 'thap',"IgnoreCase",true))
            roitab.Perturbation{i} = 'thapsigargin';
        end
        
        % Dose
        if(~(roitab.Perturbation{i} == "control"))
            dose = str2double(regexp(roitab.Folderpath{i}, '([\d.]+)(?=\uM)', 'match'));
    
            if(isempty(dose))
                dose = str2double(regexp(roitab.Folderpath{i}, '([\d.]+)(?=\µM)', 'match'));
            end
            roitab.Dose(i) = dose;
        end

    catch
        invalid_rows = [invalid_rows; i];
    end
end

roitab(invalid_rows,:) = [];

%% Iterate over the tif-stacks and count number of cells in green and yellow channel

% Select only rows where yellow channel is Ki76
subtab = roitab(logical((roitab.Yellow_ch == "Ki67") + (roitab.Yellow_ch == "ki67")),:);

% Create empty variables for storing counts
green_counts = zeros(height(subtab), 1);
yellow_counts = zeros(height(subtab), 1);

stacktable_const = parallel.pool.Constant(subtab);

% Start parallel pool
pool = gcp(); 

parfor i=1:height(subtab)
    % Get tif info from file
    path = stacktable_const.Value.Folderpath{i};
    tif_info = imfinfo(path);
    
    % Extract yellow channel (Ki67)
    yellow_ch = imread(path, 8, 'Info', tif_info);

    % Extract green channel (Stem121)
    green_ch = imread(path, 3, 'Info', tif_info);

    % Calculate the number of cells in each channel
    [count_green, imf, filtered_peaks, mask] = count_cells(green_ch(:,:,2));
    imshowpair(imf,filtered_peaks)
    hold on

    % Use the mask from green channel to only count yellow spots if
    % co-occuring with green spot
    [count_yellow, imf, filtered_peaks, ~] = count_cells(im2gray(yellow_ch), mask);
    imshowpair(imf,filtered_peaks)

    % Save counts
    green_counts(i) = count_green;
    yellow_counts(i) = count_yellow;
    i
end

% Add columns to table
subtab.Green_count = green_counts;
subtab.Yellow_count = yellow_counts;
subtab.Yellow_green_ratio = (yellow_counts ./ green_counts) .* 100;

%% Figure 5F
% Boxplot of Ki67+ cells vs Stem121+ cells (core and edge merged) for 3013, 3180 and 3054 control, dasatinib and thapsigargin

idx = logical((subtab.Perturbation == "control") | ...
    ((subtab.Perturbation == "thapsigargin") & (subtab.Dose == 4)) | ...
    ((subtab.Perturbation == "dasatinib") & (subtab.Dose == 32)));

subsub = subtab(idx,:);
subsub.Area = categorical(subsub.Area);
subsub.HGCC = categorical(subsub.HGCC); 

b = boxchart(subsub.HGCC, subsub.Yellow_green_ratio, 'GroupByColor', subsub.Perturbation);
b(3).SeriesIndex = 4;
legend
ax=gca;
ax.Box = 1;

fontsize('scale', 1.5)
ylabel('% Ki67+ cells/ Stem121+ cells')
title('Treatment effect on percent proliferating cells')

%% Test significance
cellines =  unique(subsub.HGCC);
p_values=[];
for i=1:length(cellines)
    hgcc = cellines(i);

    tbl = subsub(logical(subsub.HGCC == hgcc),:);
    tbl.Perturbation = categorical(tbl.Perturbation); 
    tbl.Perturbation = reordercats(tbl.Perturbation, {'control', 'dasatinib', 'thapsigargin'});

    [p,~,stats] = anova1(tbl.Yellow_green_ratio, tbl.Perturbation);
    multcompare(stats)
    p_values = [p_values; p];
end

%% Figure 5E
% Boxplot of Ki67+ cells vs Stem121+ cells core vs edge for 3013 and 3054

figure;
idx = (subtab.Perturbation == "control");

subsub = subtab(idx,:);

idx2 = logical((subsub.HGCC == "U3013MG") + (subsub.HGCC == "U3054MG"));

subsub2 = subsub(idx2,:);

subsub2.Area = categorical(subsub2.Area);
subsub2.HGCC = categorical(subsub2.HGCC); 

b = boxchart(subsub2.HGCC, subsub2.Yellow_green_ratio, 'GroupByColor', subsub2.Area);
b(1).SeriesIndex = 3;
b(2).SeriesIndex = 5;
legend
ax=gca;
ax.Box = 1;
ax.YLim = [0 65];

fontsize('scale', 1.5)
ylabel('% Ki67+ cells/ Stem121+ cells')
title('Percent proliferating cells in core and edge regions')

%% Test significance of core and edge comparison
p_values=[];
cellines = unique(subsub2.HGCC);
for i=1:length(cellines)
    subsub3 = subsub2(subsub2.HGCC == cellines(i),:);
    [p, tbl, stats] = anova1(subsub3.Yellow_green_ratio, subsub3.Area);
    p_values = [p_values p];
end

%% Boxplot of Ki67+ cells vs Stem121+ cells in 3013 (dasatinib dose response)
figure;
idx = ((subtab.Perturbation == "dasatinib") | (subtab.Perturbation == "control")) & (subtab.HGCC == "U3013MG");

subsub = subtab(idx,:);
subsub.Area = categorical(subsub.Area);
subsub.HGCC = categorical(subsub.HGCC);
subsub.Dose = categorical(subsub.Dose);

b = boxchart(subsub.Dose, subsub.Yellow_green_ratio);
ax=gca;
ax.Box = 1;
ax.YLim = [0 80];

% b.BoxFaceColor = [0.1 0.1 0.1];

fontsize('scale', 1.5)
ylabel('% Ki67+ cells/ Stem121+ cells')
title('Dasatinib dose response for U3013MG')

%% Boxplot of Ki67+ cells vs Stem121+ cells in 3013 (thapsigargin dose response)
figure;
idx = ((subtab.Perturbation == "thapsigargin") | (subtab.Perturbation == "control")) & (subtab.HGCC == "U3013MG");

subsub = subtab(idx,:);
subsub.Area = categorical(subsub.Area);
subsub.HGCC = categorical(subsub.HGCC);
subsub.Dose = categorical(subsub.Dose);

b = boxchart(subsub.Dose, subsub.Yellow_green_ratio);
ax=gca;
ax.Box = 1;
ax.YLim = [0 80];

fontsize('scale', 1.5)
ylabel('% Ki67+ cells/ Stem121+ cells')
title('Thapsigargin dose response for U3013MG')

