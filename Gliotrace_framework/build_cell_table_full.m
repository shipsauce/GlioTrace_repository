function outs = build_cell_table_full(tbl_in_full)
% This function takes a stacktable and iterates over the tracked cells in
% each ROI. For each track (cell), it calculates the turning angle distribution
% (TAD) and the average speed across the track. The information on all
% cells across all ROIs is saved in a new table where each row is a cell
% and columns are variables containing the following:
%  set - set number
%  exp - experiment number
%  roi - roi number (in experiment)
%  hgcc - celline
%  perturbation - treatment/control
%  dose - dose of treatment
%  morphology - dominating morphology label of viterbi track
%  tme - dominating tme label of raw track
%  trax - x-coordinates of track
%  tray - y-coordinates of track
%  noisy_labels - raw morphology label sequence
%  viterbi_path - HMM corrected path
%  tme_labels - tme label sequence
%  tad - turning angle distribution (a vector for each subsequence)
%  tad_tme - same as above but calculated for track indices based on 
%       dominating tme label
%  dominating morph idx - indicies of subsequences
%  average speed - average speed in track (one value for each subsequence)
%  average speed_tme - same as above but calculated based on dominating tme
%       label indices
%
% Input parameters:
%   tbl_in_full - stacktable
%
% Output parameters:
%   outs - cell array of cell tables for each celline
%
% @authors: Madeleine Skeppås
% @date: 140624

% Define varnames of table
varnames = {'set', 'exp', 'roi', 'hgcc', 'perturbation', 'dose', ...
'morphology', 'tme','trax', 'tray', 'noisy_labels', 'viterbi_path', ...
'tme_labels', 'tad', 'tad_tme','dominating_morph_idx', 'average_speed', 'average_speed_tme'};

cellines = unique(tbl_in_full.HGCC);
outs = struct;
out={};
info=[];

for p=1:length(cellines)
    tbl_out=table;
    subtab=tbl_in_full(tbl_in_full.HGCC == string(cellines{p}), :);

    perturbations = unique(subtab.perturbation);

    for ii=1:length(perturbations)
        pert_tab = subtab(subtab.perturbation == string(perturbations{ii}),:);

        doses = unique(pert_tab.dose);

        for jj=1:length(doses)
            tbl_in = pert_tab(pert_tab.dose == doses(jj),:);

            % Iterate over each row in the stacktable
            for i=1:height(tbl_in)
                row = tbl_in(i,:);
                paths = row.Viterbi_paths{:}; % Get HMM corrected label sequences
                traxs = row.traxs{:}; % Get x-coordinates
                trays = row.trays{:}; % Get y-coordinates
                props = row.props{:}{7}; % Get raw morph label sequences
                props_tme = row.props{:}{8}; % Get TME label sequences
            
                % Iterate over the Viterbi paths (one path=one cell)
                for j=1:length(paths)
                    path = paths{j};
                    x = traxs(~isnan(traxs(:,j)),j);
                    y = trays(~isnan(trays(:,j)),j);
                    morphology = mode(path); % Get dominating morphology of path
                    tme_dom = mode(props_tme(:,j)); % Get dominating TME association
            
                    % Extract subsequences of the path where the label agrees with the
                    % dominating label of the path
                    idx = extractSubsequences(path, morphology);
                    idx2 = extractSubsequences(props_tme(~isnan(props_tme(:,j)),j), tme_dom);
                    
                    % Calculate TAD and speed for each of the subsequences
                    tads = [];
                    speeds = [];
                    for m=1:length(idx)
                        index = idx{m};
                        [tad] = calculate_TAD(x(index), y(index));
                        tads = [tads tad];
                        speed = calculate_cell_speed(x(index), y(index), tbl_in.delta_t(i));
                        speeds = [speeds; speed];
                    end
            
                    % Do the same for dominating TME label
                    tads_tme = [];
                    speeds_tme = [];
                    for m=1:length(idx2)
                        index = idx2{m};
                        [tad] = calculate_TAD(x(index), y(index));
                        tads_tme = [tads_tme tad];
                        speed = calculate_cell_speed(x(index), y(index), tbl_in.delta_t(i));
                        speeds_tme = [speeds_tme; speed];
                    end
            
                    tbl_out = [tbl_out; array2table({row.set, row.exp, row.roi, row.HGCC row.perturbation{:} row.dose morphology tme_dom traxs(:,j) trays(:,j) props(:,j) path props_tme(:,j) tads' tads_tme' idx speeds speeds_tme},"VariableNames",varnames)];
                end
            end
            tbl_out.Properties.VariableNames = varnames;
            out{ii,jj,p} = tbl_out;
            info{ii,jj,p} = [cellines{p} ' ' perturbations{ii} ' ' num2str(doses(jj)) ' uM'];
            tbl_out = table;
        end
    end
    fprintf(['Building cell table for cell line: ' cellines{p} '...\n'])
end

outs.cell_statistics = out;
outs.info = info;

end
