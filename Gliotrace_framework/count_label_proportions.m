function [tbl_ext, data, legends] = count_label_proportions(tbl)
% This function calculates the proportions of different cell labels for
% each stack/ROI and appends them to the input table.
%
% Input parameters:
%   tbl - stacktable containing viterbi paths
%
% Output parameters:
%   tbl_ext - table extended with statistics on dominating morphologies
%   
% @authors: Madeleine Skeppås
% @date: 140624
    
majority_morphs = [];
majority_morphs_tme = [];

% Iterate through the stack. For each ROI, calculate the dominating
% identity/TME label across the time-series for each cell.
for i=1:height(tbl)
    props = tbl.props{i}{8}; % Retrieve the sequences
    idx = sum(~isnan(props),1) > 3; % Filter out sequences shorter than 4 frames
    props = props(:,idx);
    modeResults = mode(props,1); % Calculate dominating label across time

    % Use histcounts to calculate the prevalence of different labels
    varnames = {'Microglia colocalized','Vessel associated', 'Non-associated'};
    stats = array2table(histcounts(modeResults, 1:4), 'VariableNames',varnames);
    majority_morphs_tme = [majority_morphs_tme; stats];
end

% Now do the same for morphology labels
for i=1:height(tbl)
    % Same as above but using HMM corrected morphology labels
    viterbi_paths = tbl.Viterbi_paths{i};

    % Filter out sequences shorter than 4 frames
    viterbi_paths = viterbi_paths(cell2mat(cellfun(@length, viterbi_paths, 'UniformOutput', false)) > 3);

    modeResults = cell2mat(cellfun(@mode, viterbi_paths, 'UniformOutput', false)); % Get dominating label across time
    varnames = {'Branching','Diffuse translocation','Junk','Locomotion',['Perivascular' ' translocation'],'Round'};
    stats = array2table(histcounts(modeResults, 1:7), 'VariableNames',varnames); % Count the prevalence of cells with different labels
    majority_morphs = [majority_morphs; stats];
end

% Add information on the prevalence of different labels to stacktable
tbl_ext = [tbl majority_morphs majority_morphs_tme];

morphes = {};
tme_labs={};
dose_count = [];
dose_list = {};

perturbation = unique(tbl_ext.perturbation);

% For each perturbation, loop through every combination of
% dose and celline and calculate the proportion of cells in each category
% as a mean across all ROIs for that combination.
fprintf(['Looping through every combination of perturbation, dose and cell' ...
    ' line and calculating proportions of cells in each category...\n'])
for j=1:length(perturbation)
    pert = perturbation{j};
    tbl_ext_pert = tbl_ext(tbl_ext.perturbation == string(pert),:);
    dosez = unique(tbl_ext_pert.dose);
    dose_count(j) = length(dosez);
    dose_list{j} = dosez;

    % Loop through doses
    for k=1:length(dosez)
        dose_curr = dosez(k);
        tbl_ext_pert_dose = tbl_ext_pert(tbl_ext_pert.dose == dose_curr,:);
        cellines = unique(tbl_ext_pert_dose.HGCC);
    
        % Loop through each cell line
        for i=1:length(cellines)
            hgcc = cellines{i};
            tab = tbl_ext_pert_dose(strcmp(tbl_ext_pert_dose.HGCC,hgcc),:);
            
            % Calculate the ratio of different labels for each ROI
            % First look at TME labels
            try % Column indices assume HMM code has been run...
                cellcount = table2array(sum(tab(:,27:29),2));
                lab_counts = table2array((tab(:,27:29))) ./ cellcount;
            catch %...otherwise try this
                try
                    cellcount = table2array(sum(tab(:,23:25),2));
                    lab_counts = table2array((tab(:,23:25))) ./ cellcount;
                catch
                    cellcount = table2array(sum(tab(:,26:28),2));
                    lab_counts = table2array((tab(:,26:28))) ./ cellcount;
                end
            end

            % Save the mean ratios across all ROIs
            tme_labs = [tme_labs; nanmean(lab_counts,1)];
            
            % Now look at morphological labels
            try
                cellcount = table2array(sum(tab(:,[27 28 30 31 32]),2));
                lab_counts = table2array((tab(:,[27 28 30 31 32]))) ./ cellcount;
            catch
                cellcount = table2array(sum(tab(:,[26 27 29 30 31]),2));
                lab_counts = table2array((tab(:,[26 27 29 30 31]))) ./ cellcount;
            end
    
            % Save the mean ratios across all ROIs
            morphes = [morphes; nanmean(lab_counts,1)];
       
            % Save information about how many mice (exp) and ROIs were used to
            % calculate the ratios of this combination of pert, dose, celline
            legendz{i} = sprintf(['Cell line: ' hgcc ' Perturbation:  ' pert ' Dose:  ' num2str(dose_curr) ' (n = ' num2str(length(unique(tab.exp))) ')' ' (ROIs = ' num2str(height(tab)) ')']);
        end
    
        data{j,k,1} = morphes;
        data(j,k,2) = tme_labs;
        legends{j,k} = legendz;
        morphes = {};
        tme_labs = {};
        legendz = {};
    end

end
end