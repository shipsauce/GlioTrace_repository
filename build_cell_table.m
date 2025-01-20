function tbl_out = build_cell_table(tbl_in)
% This function takes a stacktable and iterates over the tracked cells in
% each ROI. For each track (cell), it calculates the turning angle distribution
% (TAD) and the average speed across the track. The information on all
% cells across all ROIs is saved in a new table where each row is a cell
% and columns are variables containing the following:
% set - set number
% exp - experiment number
% roi - roi number (in experiment)
% hgcc - celline
% perturbation - treatment/control
% dose - dose of treatment
% morphology - dominating morphology label of viterbi track
% trax - x-coordinates of track
% tray - y-coordinates of track
% noisy_labels - raw morphology label sequence
% viterbi_path - HMM corrected path
% tme_labels - tme label sequence
% tad - turning angle distribution (a vector for each subsequence)
% dominating morph idx - indicies of subsequences
% average speed - average speed in track (one value for each subsequence)
%
% Input parameters:
%   tbl_in - stacktable
%
% Output parameters:
%   tbl_out - cell table
%
% @authors: Madeleine Skeppås
% @date: 140624

% Define varnames of table
varnames = {'set', 'exp', 'roi', 'hgcc', 'perturbation', 'dose', ...
'morphology', 'trax', 'tray', 'noisy_labels', 'viterbi_path', ...
'tme_labels', 'tad', 'dominating_morph_idx', 'average_speed'};

tbl_out=table;

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
        morphology = mode(path); % Get dominating morph of path

        % Extract subsequences of the path where the label agrees with the
        % dominating label of the path
        idx = extractSubsequences(path, morphology);
        
        % Calculate TAD and speed for each of the subsequences
        tads = [];
        speeds = [];
        for m=1:length(idx)
            index = idx{m};
            [tad] = calculate_TAD(x(index), y(index));
            tads = [tads tad];
            speed = calculate_cell_speed(x(index), y(index), tbl_in.delta_t(m));
            speeds = [speeds; speed];
        end

        tbl_out = [tbl_out; array2table({row.set, row.exp, row.roi, row.HGCC row.perturbation{:} row.dose morphology traxs(:,j) trays(:,j) props(:,j) path props_tme(:,j) tads' idx speeds},"VariableNames",varnames)];
    end
    i
end

tbl_out.Properties.VariableNames = varnames;

end