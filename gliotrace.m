function gliotrace_output = gliotrace(stackfile, output_path)
% The GlioTrace framework takes a set of stabilized RGB image stacks and 
% performs cell identification, tracking and classification. The resulting
% cell tracks and labels are used to fit parameters of a Hidden Markov
% model to estimate cell morphological transitions across different cell
% lines. Based on label sequences, the proportions of cells in different states are
% calculated for each stack, and additional single cell statistics are created to
% further map the differences in migration between different morphological
% states. 
%
% Further additional statistics are average pixel movement and increase in GFP as a
% proxy for cell proliferation, measured at stack-level.
%
% A combination of single-cell statistics and stack-level statistics are
% used to estimate dose-dependent drug effects by fitting an LME.
% 
% @authors: Madeleine Skeppås
% @date: 01032025

% Perform cell tracking and classification, calculate ROI-level statistics
% and generate corresponding videos saved into the defined output path (if
% provided)
if(nargin>1)
    [slice_statistics, vasculature_statistics] = build_statistics_v3(stackfile, output_path);
else
    [slice_statistics, vasculature_statistics] = build_statistics_v3(stackfile);
end

% Force control samples to have zero dose
slice_statistics.dose(slice_statistics.perturbation == "control") = 0;

% Remove ROIs without tracked objects
idx = cellfun(@(x) isequal(x, []), slice_statistics.traxs);
tbl = slice_statistics(~idx,:);

% Fit and apply a Hidden Markov Model to cell tracks
tbl_fit = fit_apply_hmm_v3(tbl);

% Add statistics on number of cells in each class
slice_statistics_append = count_label_proportions(tbl_fit);

% Create table of single cell statistics
cell_table_full = build_cell_table_full(slice_statistics_append);

% Create table of statistics on drug effects
drug_statistics = compute_drug_effects(slice_statistics_append);

gliotrace_output = struct;
gliotrace_output.slice_statistics = slice_statistics_append;
gliotrace_output.vasculature_statistics = vasculature_statistics;
gliotrace_output.cell_statistics = cell_table_full;
gliotrace_output.drug_statistics = drug_statistics;

fprintf('Done!\n')

end