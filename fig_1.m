%% Fig 1
% This script produces material seen in figure 1 of the paper.
%
% @authors: Madeleine Skeppås
% @date: 15012025
%% Quantify the growth
% Fetch metadata about experiments
metadata=readtable('/Volumes/MyGroups$/Iron/konfokalmikroskop/Hitesh Montage and Overlays/hitesh_metadata.xlsx');

% Build table of ROIs from file
stacktable=build_stack_table(metadata,'growth_stacks_complete.txt');

% Select the relevant perturbations & dose
subtable = stacktable(logical((stacktable.perturbation == "control") + (stacktable.perturbation == "dasatinib") + (stacktable.perturbation == "thapsigargin")),:);
dose_idx = (subtable.dose == 0) + (subtable.dose == 32) + (subtable.dose == 4);
subtable = subtable(logical(dose_idx),:);
dasat_idx = (subtable.dose == 4) & (subtable.perturbation == "dasatinib");
subtable = subtable(~dasat_idx,:);

% Create empty variable for storing total intensity of green channel in 
% each frame across stack
sum_green = {};


%% Iterate through the stacks
for i=1:height(subtable)
    % Load the stack
    stack = load(subtable.file{i}); 
    stack = stack.stack;
    
    % Save green and red channel in separate variables
    gbm=stack.Tstack;
    vasc=stack.Vstack;

    % Sum over first two dimensions of stack (x,y) and save the
    % time-resolved measure of intensity (green channel) to variable
    sum_green{i} = squeeze(mean(mean(gbm)));
    i
end

% Append intensity measures to stack table
subtable.growth = sum_green';

%% Save table
growth_statistics = subtable;
save(['growth_statistics_' char(datetime('today')) '.mat'], 'growth_statistics')

%% Fig 1F, G, H

% ---------------------- Plot the growth ---------------------------------
% Load table
clear all
load("growth_statistics_22-Nov-2024.mat");

 % Retrieve the names of HGCC cellines
cellines = unique(growth_statistics.HGCC);

% Create color variable and legend array
col = [
  [13/255, 59/255, 102/255];    
  [61/255, 165/255, 217/255];   
  [115/255, 191/255, 184/255];  
  [185/255, 195/255, 93/255];   
  [254/255, 198/255, 1/255];
  [249/255, 87/255, 56/255]
];


leg1=[];

% Only include the ROIs which are labelled control
stacktable_ctrl = growth_statistics(growth_statistics.perturbation == "control",:);
stacktable_ctrl.interpolated_growth = cell(height(stacktable_ctrl),1);

% Find the smallest shortest imaging interval and use this to interpolate
% values according to this
deltat_min = min(stacktable_ctrl.delta_t);

% Loop through the cell lines
for i=1:length(cellines)
    hold on
    hgcc = cellines{i};
    interpolated_values = {};

    % Tab is a subset of stacktable only containing stacks from the same
    % cell line
    tab = stacktable_ctrl(strcmp(stacktable_ctrl.HGCC,hgcc),:);

    % Handle imaging gap in set 45 (exp 214, 215)
    if(hgcc == "U3054MG")
        % Take out correct set
        subtab = tab([1 2],:);
        
        % Exp 214
        old_vals = subtab.growth{1};
        interp_15_16 = interp1([1 11], [old_vals(15) old_vals(16)], 1:11);
        interp_27_28 = interp1([1 40], [old_vals(27) old_vals(28)], 1:40);
        new_vals = [old_vals(1:14)' interp_15_16 old_vals(17:26)' interp_27_28 old_vals(29:end)']';
        tab.growth(1) = {new_vals};

         % Exp 215
        old_vals = subtab.growth{2};
        interp_15_16 = interp1([1 11], [old_vals(15) old_vals(16)], 1:11);
        interp_26_27 = interp1([1 40], [old_vals(26) old_vals(27)], 1:40);
        new_vals = [old_vals(1:14)' interp_15_16 old_vals(17:25)' interp_26_27 old_vals(28:end)']';
        tab.growth(2) = {new_vals};
    end
    
    % Loop through stacks to interpolate values according to the smallest
    % deltat
    for k=1:height(tab)
        t2 = linspace(0,tab.delta_t(k)*(length(tab.growth{k})-1), length(tab.growth{k}));
        t1 = 0:deltat_min:t2(end);
        tab.interpolated_growth(k) = {interp1(t2,tab.growth{k},t1)};
    end

    % Create empty variable for storing values of growth
    values = nan(height(tab),max(cellfun(@length, tab.interpolated_growth)));

    % Loop through the stacks of the current cell line, extract 
    % and convert the values from the table
    for m=1:height(tab)
        vals = tab.interpolated_growth{m}';
        vals = vals./vals(1) * 100; % Divide by starting intensity and convert to procentual change
        values(m,1:length(tab.interpolated_growth{m})) = vals; % Save to variable
    end
    
    % Cut if the time-resolved measures are longer than 70 frames
    % if(size(values,2)>70)
    %     values = values(:,1:70);
    % end

    t=0:deltat_min:(deltat_min*(width(values)-1));

    if(hgcc == "U3028MG")
        t=t(1:71);
        values = values(:,1:71);
    end
    
    % Plot curves as mean + standard deviation
    subplot(3,1,1)
    p = stdshade(t,values, 0.1, col(i,:), [hgcc ' (n = ' num2str(length(unique(tab.exp))) ')' ' (ROIs = ' num2str(height(tab)) ')'], '-', 'none',8);
    leg1 = [leg1 p];
end

% Add figure details and adjust axes
title('Proliferation')
xlabel('t (h)')
ylabel('% GFP/baseline')
al = legend(leg1,"Location","bestoutside");
a=axis();
a(2)=120;
axis(a)
ax = gca;

ax.YScale = 'log';
ax.LineWidth = 0.7;

% -------------------------- Plot the adMAD ------------------------

% Load full statistics table
tbl_full = load("slice_statistics_03-Sep-2024.mat");
tbl_full = tbl_full.slice_statistics;

% Recreate empty line legend array
leg1=[];

% Only include the ROIs which are labelled control
tbl = tbl_full(tbl_full.perturbation == "control",:);
tbl.interpolated_admad = cell(height(tbl),1);

deltat_min = min(tbl.delta_t);

% Loop through the cell lines
for i=1:length(cellines)
    hold on
    hgcc = cellines{i};

    % Create the subset of stacks that belong to the current cell line
    tab = tbl(strcmp(tbl.HGCC,hgcc),:);

    % Handle imaging gap in set 45 (exp 214, 215)
    if(hgcc == "U3054MG")
        
        % Exp 214
        subtab = tab(tab.exp  == 214,:);
        for j=1:height(subtab)
            old_vals = subtab.adMAD{j};

            windowSize = 3;
            boxKernel = ones(1, windowSize) / windowSize;  % 1D box filter, normalized

            % Apply the box filter using convolution
            smoothedData = conv(old_vals, boxKernel, 'same'); % 'same' to maintain data length

            new_vals = interp1([1:15 25:36 75:93], smoothedData, 1:93);
            new_admad(j) = {new_vals'};
        end

        tab.adMAD(tab.exp == 214) = new_admad;

        % Exp 215
        subtab = tab(tab.exp  == 215,:);
        for j=1:height(subtab)
            old_vals = subtab.adMAD{j};

            windowSize = 3;
            boxKernel = ones(1, windowSize) / windowSize;  % 1D box filter, normalized

            % Apply the box filter using convolution
            smoothedData = conv(old_vals, boxKernel, 'same'); % 'same' to maintain data length

            new_vals = interp1([1:15 25:35 74:92], smoothedData, 1:92);
            new_admad(j) = {new_vals};
        end

        tab.adMAD(tab.exp == 215) = new_admad;
    end
    
    % Loop through stacks to interpolate values according to the smallest
    % deltat
    % figure;
    for k=1:height(tab)
        t2 = linspace(0,tab.delta_t(k)*(length(tab.adMAD{k})-1), length(tab.adMAD{k}));
        t1 = 0:deltat_min:t2(end);
        tab.interpolated_admad(k) = {interp1(t2,tab.adMAD{k},t1)};
    end

    % Create empty variable for storing adMAD values
    values = nan(height(tab),max(cellfun(@length, tab.interpolated_admad)));

    % Loop through the stacks of the current cell line, extract 
    % the values and remove outliers
    for m=1:height(tab)
        vals = tab.interpolated_admad{m}';
        [~, tfrm] = rmoutliers(vals, 'quartiles'); % Identify outliers
        vals(tfrm) = mean(vals); % Replace outliers with mean
        values(m,1:length(tab.interpolated_admad{m})) = vals*100;
    end
    
     % Cut if the time-resolved measures are longer than 70 frames
    % if(size(values,2)>70)
    %     values = values(:,1:70);
    % end

    t=0:deltat_min:(deltat_min*(width(values)-1));

    % Plot curves as mean + standard deviation
    subplot(3,1,2)
    hold on
    p = stdshade(t,values, 0.1, col(i,:), [hgcc ' (n = ' num2str(length(unique(tab.exp))) ')' ' (ROIs = ' num2str(height(tab)) ')'], '-', 'none',8);
    leg1 = [leg1 p];
end

% Add figure details and adjust axes
title('Adjusted mean of abs. differences (adMAD)')
xlabel('Δt (h)')
a=axis();
a(2)=90;
axis(a)
ax = gca;

% Set the axis label and tick mark colors to white
ax.XColor = [0 0 0];  
ax.YColor = [0 0 0];  

bl = legend(leg1, 'Location', 'bestoutside');
box(ax,'on')
ax.LineWidth = 0.7;

% ------------------------- Calculate and plot MSD -----------------------
% Use the control slice statistics table

% Recreate empty line legend array
leg1=[];

% Only include the ROIs which are labelled control
tbl = tbl_full(tbl_full.perturbation == "control",:);


% Loop through the cell lines
for i=1:length(cellines)
    hgcc = cellines{i};
    
    % Create the subset of stacks that belong to the current cell line
    tab = tbl(tbl.HGCC == string(hgcc),:);
    max_frames = max(tab.frames);
    deltat_max = max(tab.delta_t);

    t_max = 0:deltat_min:(max_frames*deltat_max);

    setz = unique(tab.set);

    sd_celline=nan(length(t_max),sum(cellfun(@(x) size(x,2), tab.traxs)));

    for k=1:(length(setz))
        zet = setz(k);
        subtable = tab(tab.set == zet,:);

        delta_t = subtable.delta_t(1);

        max_interpolated_length = length(0:deltat_min:(delta_t*subtable.frames(1)));
        
        % Create empty arrays to store all tracks from all stacks, based on the
        % maximum number of frames in a stack, and the total number of tracks
        % from all stacks. Arrays follow the same structure as structures traxs
        % and trays, but are meant to store all such structures together.
        % Rows: frames
        % Columns : tracks
        all_traxs = nan(max(subtable.frames), sum(cellfun(@(x) size(x,2), subtable.traxs)));
        all_trays = nan(max(subtable.frames), sum(cellfun(@(x) size(x,2), subtable.traxs)));
        
        % Handle the first stack separately
        xs = subtable.traxs{1};
        ys = subtable.trays{1};
    
        % Append the first tracks
        all_traxs(1:height(xs),1:width(xs)) = xs;
        all_trays(1:height(xs),1:width(xs)) = ys;
    
        % This is where the appended tracks end
        tailend = width(xs);
                
        % Iterate over the remaining stacks, appending tracks in the same way
        for kk=2:height(subtable)
            xs = subtable.traxs{kk};
            ys = subtable.trays{kk};
    
            all_traxs(1:height(xs),tailend+1:tailend+width(xs)) = xs;
            all_trays(1:height(xs),tailend+1:tailend+width(xs)) = ys;
    
            tailend = tailend + width(xs);
        end

        % Create empty array for squared displacement of tracks
        sd = nan(size(all_traxs));
        
        % Iterate over all tracks (column by column)
        for m=1:width(all_traxs)
            reference_pos_x=all_traxs(min(find(~isnan(all_traxs(:,m)))),m); % Find starting pos x of track
            reference_pos_y=all_trays(min(find(~isnan(all_trays(:,m)))),m); % Find starting pos y of track
            sd(1:length(find(~isnan(all_traxs(:,m)))),m) = abs(all_traxs(find(~isnan(all_traxs(:,m))),m) - reference_pos_x).^2 + ...
                abs(all_trays(find(~isnan(all_trays(:,m))),m) - reference_pos_y).^2; % Calculate how far from the starting point each coordinate in track is
        end

        % Remove tracks shorter than some threhsold
        idx = sum(~isnan(sd),1) <= 5;
        sd = sd(:, ~idx);

        sd_interpolated_set = nan(max_interpolated_length, width(sd));

        % Loop through sd values to interpolate according to the smallest deltat
        for m=1:width(sd)
            track = sd(~isnan(sd(:,m)),m);
            tracklength = sum(~isnan(track));
            t2 = linspace(0,delta_t*(tracklength-1), tracklength);
            t1 = 0:deltat_min:t2(end);
            try
                interpolated_values = interp1(t2,track,t1);
                sd_interpolated_set(1:length(interpolated_values),m) = interpolated_values;
            catch
                sd_interpolated_set(1:tracklength,m) = track;
            end
        end

        % Append to the cell line table collecting all squared displacement
        % values
        start_idx = min(find(isnan(nanmean(sd_celline))));
        sd_celline(1:height(sd_interpolated_set), start_idx:start_idx + width(sd_interpolated_set)-1) = sd_interpolated_set;
    end

    % Convert squared displacement to um
    sd_um = sd_celline .* 1.85;

    t=0:deltat_min:(deltat_min*(height(sd_um)-1));
  
    % Plot mean + standard deviation (bold line thus becomes MSD)
    subplot(3,1,3)
    hold on
    p = stdshade(t,sd_um', 0.1, col(i,:), [hgcc ' (n = ' num2str(length(unique(tab.exp))) ')' ' (ROIs = ' num2str(height(tab)) ')'], '-','none',8);
    leg1 = [leg1 p];
    ax=gca;
    ax.LineWidth = 0.7;

    i
end

% Add figure details and adjust axes
cl = legend(leg1, 'Location', 'bestoutside');
title('Cell migration (MSD)')
xlabel('t (h)')
ylabel('MSD (µm)^2')
fontsize('scale',1.5)
a=axis();
a(2)=56;
axis(a)

% Normal diffusion
D = 5;
alpha = 1;
msd_norm = 4 * D * t .^ alpha;

plot(t,msd_norm, 'LineStyle',':', 'Color','k', 'LineWidth', 2)

set(gca, 'XScale', 'log', 'YScale', 'log')

%% Supplementary figure S1

% Build TIF table from textfile
tab=readtable('fig1_viability_stainings.txt','ReadVariableNames',false,'delimiter','\t');
tab.Properties.VariableNames{1} = 'Folderpath';
tab.Day = zeros(height(tab),1);
tab.Channel = cell(height(tab),1);
tab.Well = cell(height(tab),1);

% Extract metadata on files

invalid_rows = [];

for i=1:height(tab)

    try
        % Day
        tab.Day(i) = str2double(extractBetween(tab.Folderpath{i}, 'Day_', '/'));

        % Channel
        tab.Channel{i} = regexprep(char(regexp(tab.Folderpath{i}, '/([^/]+)(?=/[^/]*$)', 'match')), '/', '');
        
        % Well
        tab.Well{i} = char(regexp(tab.Folderpath{i}, '([a-zA-Z]?\d+)(?=\.tif)', 'match'));
        
        if(length(char(regexp(tab.Folderpath{i}, '([a-zA-Z]?\d+)(?=\.tif)', 'match'))) == 1)
            tab.Well{i} = [];
        end

        if(any(cellfun(@isempty, table2cell(tab(i, :)))))
            invalid_rows = [invalid_rows; i];
        elseif(contains(tab.Folderpath{i},'results', 'IgnoreCase',true))
            invalid_rows = [invalid_rows; i];
        end

    catch
        invalid_rows = [invalid_rows; i];
    end

end

tab(invalid_rows,:) = [];
tab((tab.Channel == "Green" | tab.Channel == "Red-Green-Blue"),:) = [];

Cell_count = zeros(height(tab), 1);

stacktable_const = parallel.pool.Constant(tab);

% Start parallel pool
pool = gcp(); 

% Calculate the number of cells in blue and red channels
parfor i=1:height(tab)
    % Get path
    path = stacktable_const.Value.Folderpath{i};
    
    % Read TIF
    im = imread(path);

    if(stacktable_const.Value.Channel{i} == "Blue")
        im=im(:,:,3);
    else
        im=im(:,:,1);
    end

    % Calculate the number of cells
    [count, imf, filtered_peaks, ~] = count_cells(im);

    % Save count
    Cell_count(i) = count;
    i
end

tab.Cell_count = Cell_count;

% Calculate % of Ethd1+/Dapi+ cells

tab_sorted = sortrows(tab, ["Day", "Well"]);

for i=1:height(tab_sorted)
    if(mod(i,2))
        % For odd rows (blue channel)
        tab_sorted.Ratio_dead_cells(i) = (tab_sorted.Cell_count(i) / tab_sorted.Cell_count(i)) * 100; % 1 
    else
        % For even rows (red channel)
        tab_sorted.Ratio_dead_cells(i) = (tab_sorted.Cell_count(i) / tab_sorted.Cell_count(i-1)) * 100; % Cell count red channel/ cell count blue channel
    end
end

% Bar chart

tbl_red = tab_sorted(tab_sorted.Channel == "Red",:);
tbl_red(tbl_red.Day == 8,:) = [];
tbl_red.Day = categorical(tbl_red.Day);

col = [
  [115/255, 191/255, 184/255];  % First color (Soft Teal)
  [70/255, 130/255, 180/255];   % Second color (Steel Blue)
  [255/255, 99/255, 71/255];    % Third color (Tomato Red)
  [189/255, 183/255, 107/255];  % Fourth color (Olive Drab)
];

% Plot violin plot of percentage of cells interacting with microglia
figure;
dabarplot(tbl_red.Ratio_dead_cells,'groups', tbl_red.Day,'colors', col);
ax=gca;
ax.YLim = [0 100];
ax.Box = 1;

fontsize('scale', 1.5)

ylabel('% Ethd1+/Dapi+ cells')
xlabel('Day')
title('Percentage of dead cells')


