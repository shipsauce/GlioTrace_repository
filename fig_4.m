% Fig 4
% This script produces all the material seen in figure 4 of the paper.

%% Add calculations of speed and average adMAD to the table

% Generate table with additional information on morphology ratio
tbl_fit_orig = tbl_fit;
mode = "mean";
morph = "morph";
perturbation = {"control"};
varnames = {'Branching','Diffuse translocation','Junk','Locomotion',['Perivascular' ' translocation'],'Round'};
tbl_fit = plot_proportions(tbl_fit, mode, perturbation, varnames,morph);


% Add empty columns to table
tbl_fit.alpha = cell(height(tbl_fit),1);
tbl_fit.D = cell(height(tbl_fit),1);
tbl_fit.speed = cell(height(tbl_fit),1);
tbl_fit.adMAD_mean = cell(height(tbl_fit),1);
tbl_fit.prolif_mean = cell(height(tbl_fit),1);


% Force control samples to have zero dose
tbl_fit.dose(tbl_fit.perturbation == "control") = 0;

% Define the cellines
cellines = unique(tbl.HGCC);

for i=1:height(tbl_fit)
    % Add speed
    tbl_fit.speed{i} = calculate_cell_speed(tbl_fit.traxs{i}, tbl_fit.trays{i},tbl_fit.delta_t(i));

    % Add avg adMAD
    tbl_fit.adMAD_mean{i} = mean(tbl_fit.adMAD{i});

    % Add avg proliferation
    tbl_fit.prolif_mean{i} = mean(tbl_fit.sum_green{i});
end

%% Calculate MSD curve parameters for every combination of cell line, perturbation & dose, set

% Iterate over the cellines
for i=1:length(cellines)
    hgcc = cellines{i};
    subtable = tbl_fit(tbl_fit.HGCC == string(hgcc),:);

    perturbations = unique(subtable.perturbation);
    
    % Iterate over the different perturbations tested in celline
    for j=1:length(perturbations)
        pert = perturbations{j};
        
        subtable_1 = subtable(subtable.perturbation == string(pert),:);

        doses = unique(subtable_1.dose);
        
        % Iterate over the doses of current perturbation
        for jj=1:length(doses)
            
            dose = doses(jj);
            
            subtable_2 = subtable_1(subtable_1.dose == dose,:);
            
            idx = logical((tbl.HGCC == string(hgcc)) .* (tbl.perturbation == string(pert)) .* (tbl.dose == dose));
            
            max_frames = max(subtable_2.frames);
            deltat_max = max(subtable_2.delta_t);
            deltat_min = min(subtable_2.delta_t);
            
            t_max = 0:deltat_min:(max_frames*deltat_max);
            
            setz = unique(subtable_2.set);
            
            sd_subset=nan(length(t_max),sum(cellfun(@(x) size(x,2), subtable_2.traxs)));
            
            % Iterate over all the sets and interpolate the values from
            % each set separately, since each set can have distinct delta t
            
            for k=1:(length(setz))
                zet = setz(k);
                subtable_3 = subtable_2(subtable_2.set == zet,:);
                
                delta_t = subtable_3.delta_t(1);
                
                max_interpolated_length = length(0:deltat_min:(delta_t*subtable_3.frames(1)));
                
                % Create empty arrays to store all tracks from all stacks, based on the
                % maximum number of frames in a stack, and the total number of tracks
                % from all stacks. Arrays follow the same structure as structures traxs
                % and trays, but are meant to store all such structures together.
                % Rows: frames
                % Columns : tracks
                all_traxs = nan(max(subtable_3.frames), sum(cellfun(@(x) size(x,2), subtable_3.traxs)));
                all_trays = nan(max(subtable_3.frames), sum(cellfun(@(x) size(x,2), subtable_3.traxs)));
                
                % Handle the first stack separately
                xs = subtable_3.traxs{1};
                ys = subtable_3.trays{1};
                
                % Append the first tracks
                all_traxs(1:height(xs),1:width(xs)) = xs;
                all_trays(1:height(xs),1:width(xs)) = ys;
                
                % This is where the appended tracks end
                tailend = width(xs);
                
                % Iterate over the remaining stacks, appending tracks in the same way
                for kk=2:height(subtable_3)
                    xs = subtable_3.traxs{kk};
                    ys = subtable_3.trays{kk};
                    
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
                start_idx = min(find(isnan(nanmean(sd_subset))));
                sd_subset(1:height(sd_interpolated_set), start_idx:start_idx + width(sd_interpolated_set)-1) = sd_interpolated_set;
            end
            
            % Convert squared displacement to um
            sd_um = sd_subset .* 1.85;

            t=0:deltat_min:(deltat_min*(height(sd_um)-1));

            % Calculate MSD
            msd = nanmean(sd_um,2);
         
            % Fit a straight line to get the diffusion 
            % parameters for the first time points
            log_y = log(msd(3:15));
            log_x = log(t([3:15]));
            p = polyfit(log_x, log_y, 1);
            log_y_fit = polyval(p, log_x);

            % Plot the fit, the MSD curve and the points used to fit the
            % line
            plot(log_x, log_y_fit, 'k-')
            hold on
            plot(log(t(3:end)),log(msd(3:end)), 'b:')
            hold on
            plot(log_x, log_y, 'ro')
            hold off
            
            % Add information on D and alpha to the relevant ROIs for
            % current combination of celline, perturbation and dose
            tbl_fit.alpha(idx) = {p(1)};
            tbl_fit.D(idx) = {p(2)};
        
        end
        
    end

    i

end
%% Fig 4B

% Select the relevant ROIs
idx = logical((tbl_fit.perturbation == "control") | ...
    ((tbl_fit.perturbation == "thapsigargin") & (tbl_fit.dose == 4)) | ...
    ((tbl_fit.perturbation == "dasatinib") & (tbl_fit.dose == 32)));
subtable = tbl_fit(idx,:);
subtable = subtable(~(subtable.HGCC == "U3028MG"),:);
subtable = sortrows(subtable,"perturbation","ascend");
subtable.HGCC = categorical(subtable.HGCC);
subtable.speed = cell2mat(subtable.speed);

% Plot boxchart of the speed of cells grouped by perturbation
b=boxchart(subtable.HGCC, subtable.speed, 'GroupByColor',subtable.perturbation);
legend
xlabel('HGCC Cell line')
ylabel('Cell speed (um/h)')
title("Effect of treatment")
fontsize('scale',1.3)
b(3).SeriesIndex = 4;

cellines = unique(subtable.HGCC);

% Calculate significance
p_values=[];
for i=1:length(cellines)
    hgcc = cellines(i);

    tbl = subtable(logical(subtable.HGCC == hgcc),:);

    [p,~,stats] = anova1(tbl.speed, tbl.perturbation);
    multcompare(stats)
    p_values = [p_values; p];
end

%% Fig 4C

% Define the relations between dependent and independent variables to be
% tested
formulas = {'speed ~ perturbation:dose + (1|set)', 'D ~ perturbation:dose + (1|set)', 'alpha ~ perturbation:dose + (1|set)', ...
    'growth_rate ~ perturbation:dose + (1|set)', 'adMAD_mean ~ perturbation:dose + (1|set)'};

no_of_stats = 5; 
size_i = length(unique(tbl_fit.HGCC))-1;
size_j = length(formulas);
size_k = no_of_stats;
stats_dasat = nan(size_i, size_j, size_k);
stats_thaps = nan(size_i, size_j, size_k);

cellines = unique(tbl_fit.HGCC);
cellines = cellines([1 3 4 5 6]);

% Loop over the cellines
for i=1:length(cellines)
    
    % Extract the relevant data and reformat
    tab=tbl_fit(logical(tbl_fit.HGCC == string(cellines{i})),:); 
    tab.Properties.VariableNames(31) = "Perivascular_translocation";
    tab.Properties.VariableNames(28) = "Diffuse_translocation";
    tab.alpha = cell2mat(tab.alpha);
    tab.adMAD_mean = cell2mat(tab.adMAD_mean);
    tab.prolif_mean = cell2mat(tab.prolif_mean);
    tab.D = cell2mat(tab.D);
    tab.HGCC = categorical(tab.HGCC);
    tab.perturbation = categorical(tab.perturbation);
    tab.speed = cell2mat(tab.speed);
    tab.perturbation = reordercats(tab.perturbation, string([{'control', 'thapsigargin', 'dasatinib'} sort(setdiff(unique(tab.perturbation), {'control', 'thapsigargin', 'dasatinib'}))']));

    % Iterate over the relations and fit a linear mixed effects model
    for j=1:length(formulas)
        lme = fitlme(tab, formulas{j});

        % Save information on estimate, SE, tStat, DF and pValue
        for k=1:no_of_stats
            stats_thaps(i,j,k) = lme.Coefficients{2,k+1}; 
            stats_dasat(i,j,k) = lme.Coefficients{3,k+1}; 
        end
    end

    i

end

%% Fig 4D

% Plot heatmap of the normalized estimate for each treatment
figure;
for iter=1:2
    subplot(1,3,iter)
    if(iter==1)
        M=stats_thaps(:,:,1)';
    elseif(iter==2)
        M=stats_dasat(:,:,1)';
    end
    I=[1 2 4 5]; % Exclude U3179MG
    M=M(:,I);
    M=M./max(abs(M'))';
    M = M([1 5 2 3 4],:);
    
    h=heatmap(M);
    h.YDisplayLabels = {'cell speed', 'adMAD', 'D', 'alpha', 'growth rate'};
    h.XDisplayLabels=cellines(I);
    colormap(redbluecmap)
    clim([-1 1]*1.2)
end

fontsize('scale', 1.5)

%% Fig 4E

% Plot proportions of morphological classes  before and after treatment
mode = "combined";
morph = "morph";
perturbation = {"control", "dasatinib", "thapsigargin"};
varnames = {'Branching','Diffuse translocation','Junk','Locomotion',['Perivascular' ' translocation'],'Round'};
plot_proportions(tbl_fit_orig, mode, perturbation, varnames,morph);
fontsize('scale', 1.5)

%% Fig 4F
% Plot differential transition hierarchies after treatment

% Select the correct rows
u=unique(tbl_fit.HGCC);
u = u([1 3 4 5 6]); % Remove U3028MG because it has not been treated
F=[];
for i=1:length(u)
    f1=strmatch('control',tbl_fit.perturbation);
    f2 = strmatch('dasatinib',tbl_fit.perturbation);
    f3 = strmatch('thapsigargin',tbl_fit.perturbation);
    f4=strmatch(u{i},tbl_fit.HGCC);
    f5=intersect(f1,f4); % ctrl
    f6=intersect(f2,f4); % dasatinib
    f7=intersect(f3,f4); % thapsigargin
    F(i,1)=f5(1);
    F(i,2)=f6(1);
    F(i,3)=f7(1);
end

M={};

% Extract data in A_est from HMM fitting
for i=1:height(F)
    M(:,:,1,i)=tbl_fit.A_est(F(i,1)); % ctrl
    M(:,:,2,i)=tbl_fit.A_est(F(i,2)); % dasatinib
    M(:,:,3,i)=tbl_fit.A_est(F(i,3)); % thapsigargin
end

state_names = {'B','D','Junk','L','P','R'};

% Plot Dasatinib subplot
for i=1:height(F)
    subplot(3,2,i)
    tittle = ['HGCC: ' u{i} ' Treatment: Dasatinib 32 um'];
    plot_transition_differences_v2(M(:,:,1,i), M(:,:,2,i), state_names, tittle)
end

sgtitle('Treatment effect on state transition hierarchies (TREATMENT - CTRL)')
fontsize('scale', 1.5)

figure;
% Plot Thapsigargin subplot
for i=1:height(F)
    subplot(3,2,i)
    tittle = ['HGCC: ' u{i} ' Treatment: Thapsigargin 4 um'];
    plot_transition_differences_v2(M(:,:,1,i), M(:,:,3,i), state_names, tittle)
end

sgtitle('Treatment effect on state transition hierarchies (TREATMENT - CTRL)')
fontsize('scale', 1.5)

%% Supplementary Fig S4
% Dose response thapsigargin U3013MG
idx = logical((tbl_fit.perturbation == "control") + (tbl_fit.perturbation == "thapsigargin")) .* (tbl_fit.HGCC == "U3013MG");
tab = tbl_fit(logical(idx),:);
tab.HGCC = categorical(tab.HGCC);
tab.speed = cell2mat(tab.speed);
tab.dose = categorical(tab.dose);

boxchart(tab.dose, tab.speed)
xlabel('Dose (uM)')
ylabel('Cell speed (um/h)')
title("Dose response of thapsigargin on U3013MG")
fontsize('scale',1.5)

%% Supplementary Fig S4
% Dose response dasatinib U3013MG

idx = logical((tbl_fit.perturbation == "control") + (tbl_fit.perturbation == "dasatinib")) .* (tbl_fit.HGCC == "U3013MG");
tab = tbl_fit(logical(idx),:);
tab.HGCC = categorical(tab.HGCC);
tab.speed = cell2mat(tab.speed);
tab.dose = categorical(tab.dose);

boxchart(tab.dose, tab.speed)
xlabel('Dose (uM)')
ylabel('Cell speed (um/h)')
title("Dose response of dasatinib on U3013MG")
fontsize('scale',1.5)


