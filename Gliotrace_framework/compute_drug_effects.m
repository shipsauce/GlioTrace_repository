function drug_stats = compute_drug_effects(tbl_in)
% This function takes a stacktable and adds additional statistics on cell
% speed and average adMAD and proliferation, as well as MSD parameters for
% linear fit of log-log MSD for every combination of celline, perturbation and dose.
% These statistics are used to fit a linear mixed effects model to estimate
% dose-dependent drug effects corrected for batch effects formulated as
% sets (collections of experiments). 
%
% @authors: Madeleine Skeppås
% @date: 15022025

% Add empty columns to table
tbl_in.alpha = cell(height(tbl_in),1);
tbl_in.D = cell(height(tbl_in),1);
tbl_in.speed = cell(height(tbl_in),1);
tbl_in.adMAD_mean = cell(height(tbl_in),1);
tbl_in.prolif_mean = cell(height(tbl_in),1);


% Force control samples to have zero dose
tbl_in.dose(tbl_in.perturbation == "control") = 0;

% Define the cellines
cellines = unique(tbl_in.HGCC);

for i=1:height(tbl_in)
    % Add speed
    tbl_in.speed{i} = calculate_cell_speed(tbl_in.traxs{i}, tbl_in.trays{i},tbl_in.delta_t(i));

    % Add avg adMAD
    tbl_in.adMAD_mean{i} = mean(tbl_in.adMAD{i});

    % Add avg proliferation
    tbl_in.prolif_mean{i} = mean(tbl_in.sum_green{i});
end

% Fit line to log-log MSD plot and get params
tbl_ext = fit_msd(tbl_in);

% Define the relations between dependent and independent variables to be
% tested
formulas = {'speed ~ perturbation:dose + (1|set)', 'D ~ perturbation:dose + (1|set)', 'alpha ~ perturbation:dose + (1|set)', ...
    'growth_rate ~ perturbation:dose + (1|set)', 'adMAD_mean ~ perturbation:dose + (1|set)'};

no_of_stats = 5; % Number of parameters from lme fit
no_of_perts = length(unique(tbl_ext.perturbation)) - 1; % Number of perturbations minus control
size_i = length(unique(tbl_ext.HGCC)); % Number of cell lines
size_j = length(formulas); % Number of explanatory variables to test 
size_k = no_of_stats;
size_l = no_of_perts;

stats = nan(size_j, size_k, size_l, size_i);

cellines = unique(tbl_ext.HGCC);

% Loop over the cellines
for i=1:length(cellines)
    
    % Extract the relevant data and reformat
    tab=tbl_ext(logical(tbl_ext.HGCC == string(cellines{i})),:); 
    tab.Properties.VariableNames(31) = "Perivascular_translocation";
    tab.Properties.VariableNames(28) = "Diffuse_translocation";
    tab.alpha = cell2mat(tab.alpha);
    tab.adMAD_mean = cell2mat(tab.adMAD_mean);
    tab.prolif_mean = cell2mat(tab.prolif_mean);
    tab.D = cell2mat(tab.D);
    tab.HGCC = categorical(tab.HGCC);
    tab.perturbation = categorical(tab.perturbation);
    tab.speed = cell2mat(tab.speed);
    tab.perturbation = reordercats(tab.perturbation, string([{'control'} sort(setdiff(unique(tab.perturbation), {'control'}))']));

    % Iterate over the relations and fit a linear mixed effects model
    for j=1:length(formulas)
        lme = fitlme(tab, formulas{j});

        % Save information on estimate, SE, tStat, DF and pValue
        for k=1:no_of_stats
            for l=1:no_of_perts
                try
                    stats(j,k,l,i) = lme.Coefficients{l+1,k+1}; 
                catch 
                    % For some cellines, not all perturbations have been tested
                    stats(j,k,l,i) = NaN;
                end
            end
        end
    end

    fprintf(['Fitting dose-dependent drug effects for cell line: ' cellines{i} '...\n'])

end

drug_stats = {};

for i=1:size_i % No of cellines
    for l=1:size_l % No of perturbations
        drug_stats{i,l} = array2table(stats(:,:,l,i), 'VariableNames', {'Estimate', 'SE', 'tStat', 'DF', 'pValue'}, 'RowNames',formulas);
    end
end

drug_stats = array2table(drug_stats, 'VariableNames', string(sort(setdiff(unique(tab.perturbation), {'control'}))), 'RowNames',cellines);

end