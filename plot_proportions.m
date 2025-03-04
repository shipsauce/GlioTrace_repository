function tbl_ext = plot_proportions(tbl, style, perturbation, varnames,label)
% This function generates a stacked barplot displaying the proportions
% of different cell morphologies for each celline. Dominating morphologies 
% for different cell lines (+ possible perturbations & doses) are 
% inferred using Viterbi paths.
%
% Input parameters:
%   tbl - stacktable containing viterbi paths
%   style - switch bewteen plotting mean of ROIs or all ROIs or
%               perturbations combined
%   varnames - cell array of class names
%   label - whether to plot morphology or TME label proportions
%
% Output parameters:
%   tbl_ext - table extended with statistics on dominating morphologies
%   
% @authors: Madeleine Skeppås
% @date: 140624
    
majority_morphs = [];

% Iterate through the stack. For each ROI, calculate the dominating
% identity/label across the time-series for each cell. The sequences used
% depends on label of interest.
if(label == "tme")
    for i=1:height(tbl)
        props = tbl.props{i}{8}; % Retrieve the sequences
        modeResults = mode(props,1); % Calculate dominating label across time

        % Use histcounts to calculate the prevalence of different labels
        stats = array2table(histcounts(modeResults, 1:4), 'VariableNames',varnames);
        majority_morphs = [majority_morphs; stats];
        i
    end
else
    for i=1:height(tbl)
        % Same as above but using HMM corrected morphology labels
        viterbi_paths = tbl.Viterbi_paths{i};

        % Filter out sequences shorter than 4 frames
        viterbi_paths = viterbi_paths(cell2mat(cellfun(@length, viterbi_paths, 'UniformOutput', false)) > 3);

        modeResults = cell2mat(cellfun(@mode, viterbi_paths, 'UniformOutput', false));
        stats = array2table(histcounts(modeResults, 1:7), 'VariableNames',varnames);
        majority_morphs = [majority_morphs; stats];
    end
end

% Add information on the prevalence of different labels to stacktable
tbl_ext = [tbl majority_morphs];

morphes = {};
dose_count = [];
dose_list = {};

% For each perturbation of interest, loop through every combination of
% dose and celline and 
% Loop through each perturbation
for j=1:length(perturbation)
    pert = perturbation{j};
    tbl_ext_pert = tbl_ext(tbl_ext.perturbation == pert,:);
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
            if(label == "tme")
                try % Column indices assume HMM code has been run...
                    cellcount = table2array(sum(tab(:,27:29),2));
                    lab_counts = table2array((tab(:,27:29))) ./ cellcount;
                catch %...otherwise try this
                    cellcount = table2array(sum(tab(:,23:25),2));
                    lab_counts = table2array((tab(:,23:25))) ./ cellcount;
                end
            else
                cellcount = table2array(sum(tab(:,[27 28 30 31 32]),2));
                lab_counts = table2array((tab(:,[27 28 30 31 32]))) ./ cellcount;
            end
    
            % If mean, save the mean ratios across all ROIs
            if(~(style=="all"))
                morphes = [morphes; nanmean(lab_counts,1)];
            else % Keep the information about all ROIs
                morphes{i} = lab_counts;
            end
       
            % Save information about how many mice (exp) and ROIs were used to
            % calculate the ratios of this combination of pert, dose, celline
            legendz{i} = sprintf([hgcc ' (n = ' num2str(length(unique(tab.exp))) ')' ' (ROIs = ' num2str(height(tab)) ')']);
        end
    
        data{j,k} = morphes;
        legends{j,k} = legendz;
        morphes = {};
        legendz = {};
    end

end

% Now plot
if(~(style == "combined"))
    for j=1:length(perturbation)
        
        if(~(style=="all"))
            figure
        end
    
        for k=1:dose_count(j)
            
            if(style=="mean")
        
                cols = ceil(sqrt(dose_count(j)));
                rows = ceil(dose_count(j) / cols);
                
                subplot(rows, cols , k)
        
                if(size(data{j,k},1) == 1)
                    data{j,k} = [data{j,k}{:}; zeros(size(data{j,k}{:}))];
                    bar(data{j,k}, 'stacked', 'BarWidth', 0.7, 'EdgeColor','none')
                else
                    bar(cell2mat(data{j,k}), 'stacked', 'BarWidth', 0.7, 'EdgeColor','none')
                end
        
                colors = [0.8660,    0.3290,  0;
                0.3290,    0.7130,    1.0000;
                
                0.9960,    0.5640,    0.2620;
                0.4540,    0.9210,    0.8540;
                     0,   0.6390,   0.6390
                ];
                
                ylim([0 1])
                colororder(colors);
            
                if(label == "tme")
                    legend(varnames, 'Location', 'bestoutside');
                    colororder(pink(7))
                elseif(perturbation{j} == "control")
                    legend(varnames([1 2 4 5 6]), 'Location', 'bestoutside');
                end
                xticks(1:length(legends{j,k}))
                xticklabels(legends{j,k})
                title(sprintf(['[Perturbation: ' char(perturbation{j}) ' Dose: ' num2str(dose_list{j}(k)) ' um]']))
            else
    
                colors = [0.8660,    0.3290,  0;
                0.3290,    0.7130,    1.0000;
                
                0.9960,    0.5640,    0.2620;
                0.4540,    0.9210,    0.8540;
                     0,   0.6390,   0.6390
                ];
    
                cols = ceil(sqrt(size(data{j,k},2)));
                rows = ceil(size(data{j,k},2) / cols);
                figure;
    
                tiledlayout(rows,cols,"TileSpacing","compact")
                for i=1:size(data{j,k},2)
                    nexttile(i)
                    bar(data{j,k}{i}, 'stacked', 'BarWidth', 0.7, 'EdgeColor','none')
    
                    title([legends{j,k}{i} ' [Perturbation: ' char(perturbation{j}) ' Dose: ' num2str(dose_list{j}(k)) ' um]'])
                    ylim([0 1])
                end
                sgtitle('Distribution of labels across ROIs')
                colororder(colors);
            end
        
        end
    fontsize('scale', 1.5)
    end
else
    colors = [0.8660,    0.3290,  0;
            0.3290,    0.7130,    1.0000;
            
            0.9960,    0.5640,    0.2620;
            0.4540,    0.9210,    0.8540;
                 0,   0.6390,   0.6390
            ];

            result_dasat_thaps = reshape([data{1,1}([1 3 4 5 6]), data{2,6}, data{3,4}]', [], 1);
            bar(cell2mat(result_dasat_thaps), 'stacked', 'BarWidth', 0.7, 'EdgeColor','none')

            ylim([0 1])
            colororder(colors);
            legend(varnames([1 2 4 5 6]), 'Location', 'bestoutside');
            xticks(1:15)
            xticklabels({"U3013MG ctrl", "U3013MG dasat", "U3013MG thaps", ...
                "U3054MG ctrl", "U3054MG dasat", "U3054MG thaps", ...
                "U3179MG ctrl", "U3179MG dasat", "U3179MG thaps", ...
                "U3180MG ctrl", "U3180MG dasat", "U3180MG thaps", ...
                "U3220MG ctrl", "U3220MG dasat", "U3220MG thaps"})

            ylabel("% cells")

            title("Effect of dasatinib 32um or thapsigargin 4um treatment")
            fontsize('scale', 1.5)
end
end
