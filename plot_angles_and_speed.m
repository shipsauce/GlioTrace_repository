function plot_angles_and_speed(tbl_ext, celline, perturbation)
% This function plots the turning angle distribution (TAD) and cell speed
% of individual cells from different morphological cateogories for a given
% cell line of interest. Based on information in the stacktable, a new
% table is built with build_cell_table where one row contains information
% on a single cell. This information is then separated based on the
% dominating morphology of the single cells across time and used to derive
% TAD and speed. An Anova test is calculated for the speed data.
%
% Input parameters:
%   tbl_ext - stacktable
%   celline - hgcc celline of interest
%   perturbation - perturbation of interest
%   
% @authors: Madeleine Skeppås
% @date: 140624

% Choose only the stacks which have the right combination of perturbation,
% dose and celline
tab=tbl_ext(logical((tbl_ext.HGCC == celline) .* (tbl_ext.perturbation == perturbation)),:);
tab = tab(logical((tab.dose == max(tab.dose))),:);

% Build cell table
cell_tab = build_cell_table(tab);
cell_tab.average_tad = cellfun(@mean, cell_tab.tad);

% Define the morphological categories
morph = {'Branching', 'Diffuse translocation', 'Junk','Locomotion', 'Perivascular translocation', 'Round'};

% Plot TAD
tt = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
for i=1:length(morph)
    nexttile
    tad = cell2mat(cell_tab.tad(cell2mat(cell_tab.morphology) == i,:));
    polarhistogram(tad)
    title(['Morphology: ' morph{i}])
end

% Define colors for violin plot
colors = [0.8660,    0.3290,  0;
0.3290,    0.7130,    1.0000;
0.9960,    0.5640,    0.2620;
0.4540,    0.9210,    0.8540;
     0,   0.6390,   0.6390
];

figure;

% Aggregate and log-transform data on speed
for i=1:length(morph)
    speed = cell2mat(cell_tab.average_speed(cell2mat(cell_tab.morphology) == i,:));
    data{i} = log(speed+1);
end

data = data(:,[1 2 4 5 6]);

data1=[];
group_idx = [];

% Create a grouping variable
for i=1:5
    data1 = [data1; data{i}];
    group_idx = [group_idx; repmat(i,[1 length(data{i})])'];
end

% Plot violin plot of speed
daviolinplot(data1,'groups', group_idx,'colors', colors,'outsymbol','k+',...
    'boxcolors','same','scatter',1,'jitter',1, 'scattercolors', 'same',...
    'scatteralpha', 0.4, 'violinalpha', 1, 'boxalpha', 0.9,'xtlabels', morph([1 2 4 5 6]));

% Perform ANOVA
[p, tbl, stats] = anova1(data1, group_idx);
multcompare(stats)

fontsize("scale",1.5)

end


